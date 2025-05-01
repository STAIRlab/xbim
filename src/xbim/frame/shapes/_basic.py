from shps.frame.polygon import Polygon
from shps.frame import TriangleSection
import numpy as np


class _BasicShape(Polygon):
    def __init__(self, mesh_size=None, **kwds):
        self._mesh_kwds = {}

        super().__init__(
            self.exterior(),
            self.interior(),
            mesh_size=mesh_size,
            **kwds)

    # Rendering
    def interior(self):
        return []

    # Simulation
    def elastic(self):
        return self._create_model().elastic()

    def create_fibers(self, mesh_scale=None, **kwds):
        return self._create_model(mesh_scale=mesh_scale).fibers(**kwds)

    def _create_model(self, **kwds):
        # NOTE: WideFlange overloads this to include shear warping;
        # Need to rethink how to do this generally
        mesh = self._create_mesh(**kwds)

        return TriangleSection.from_meshio(mesh,
                                           warp_shear=False)



class WideFlange(_BasicShape):
    def __init__(self, d, b, t=None, tw=None, tf=None,
                 b2=None, t2=None,
                 k = None,
                 mesh_scale=None,
                 saint_venant=None, **kwds):
        self.d  = d
        self.bf = bf = b
        if tf is None:
            tf = t
        self.tf = tf

        if b2 is None:
            b2 = b
        self.b2 = b2

        if t2 is None:
            t2 = tf
        self.t2 = t2

        if tw is None:
            tw = tf
        self.tw = tw

        self.k = k

        super().__init__(mesh_size=min(tf, tw)/3, **kwds)

    def exterior(self):
        bf = self.bf
        b2 = self.b2
        tf = self.tf
        t2 = self.t2
        tw = self.tw
        d = self.d

        y_top = d / 2.0
        y_bot = -d / 2.0

        pts = np.array([
            [-bf / 2, y_top],                 # 1  top-left flange
            [ bf / 2, y_top],                 # 2
            [ bf / 2, y_top - tf],            # 3  step down into web
            [ tw / 2, y_top - tf],            # 4
            [ tw / 2, y_bot + t2],            # 5  down web
            [ b2 / 2, y_bot + t2],            # 6  step into bottom flange
            [ b2 / 2, y_bot],                 # 7
            [-b2 / 2, y_bot],                 # 8
            [-b2 / 2, y_bot + t2],            # 9
            [-tw / 2, y_bot + t2],            # 10 up web
            [-tw / 2, y_top - tf],            # 11
            [-bf / 2, y_top - tf],            # 12
        ], dtype=float)

        return pts

#         # Area and moment of inertia
#         self.A  = tw*(d - tf - t2) + bf*tf + b2*t2

# #       Iy and Iz are wrong for tf !=  t2
# #       self.Iy = tw*(d - tf - t2)**3/12.0 + bf*tf*(0.5*(d - tf))**2 + b2*t2*(0.5*(d - t2)) ** 2
#         self.Iz = 2*tf*bf**3/12

    def shear_factor(self, nu=0.3):
        b  = self.bf
        tf = self.tf
        tw = self.tw
        d  = self.d

        m = 2*b*tf/(d*tw)
        n = b/d
        return (10*(1+nu)*(1+3*m)**2)/((12+72*m + 150*m**2 + 90*m**3) + nu*(11+66*m + 135*m**2 + 90*m**3) + 30*n**2*(m + m**2) + 5*nu*n**2*(8*m+9*m**2))


    def _create_patches(self, mesh_scale=None, mesh_kwds={}):

        from shps.frame import patch

        if mesh_scale is None:
            mesh_scale = 1/3
        bf = self.bf
        b2 = self.b2
        tf = self.tf
        t2 = self.t2
        tw = self.tw
        d  = self.d

        if self.k is None:
            quads = []
        else:
            k = self.k
            r = k - tf
            h = d / 2 - tf  # y coordinate of top of flange

            quads = [
                patch.quad(vertices=[
                    (-tw/2,     d/2- k),
                    (-tw/2,     h),
                    (-tw/2 - r, h),
                    (-tw/2 - r + r/2**0.5, d/2-k + r/2**0.5),
                ]),
                patch.quad(vertices=[
                    ( tw/2,     d/2-k),
                    ( tw/2,     h),
                    ( tw/2 + r, h),
                    ( tw/2 + r - r/2**0.5, d/2-k + r/2**0.5),
                ]),
                patch.quad(vertices=[
                    (-tw/2,    -h + r),
                    (-tw/2,    -h),
                    (-tw/2 - r, -h),
                    (-tw/2 - r + r/2**0.5, -d/2 + k - r/2**0.5),
                ]),
                patch.quad(vertices=[
                    ( tw/2,    -h + r),
                    ( tw/2,    -h),
                    ( tw/2 + r, -h),
                    ( tw/2 + r - r/2**0.5, -d/2 + k - r/2**0.5),
                ]),
            ]

        yoff = ( d - tf) / 2

        return [
            patch.rect(corners=[[-bf/2,        yoff-tf/2],[bf/2,  yoff+tf/2]]),# ,  divs=(nfl, nft), rule=int_typ),
            patch.rect(corners=[[-tw/2,       -yoff+tf/2],[tw/2,  yoff-tf/2]]),# ,  divs=(nwt, nwl), rule=int_typ),
            patch.rect(corners=[[-b2/2, -(d - t2)/2-t2/2],[bf/2, -yoff+tf/2]]),# ,  divs=(nfl, nft), rule=int_typ),
            *quads
        ]


    def _create_model(self, mesh_scale=None):
        """
        Saritas and Filippou (2009) "Frame Element for Metallic Shear-Yielding Members under Cyclic Loading"
        """
        b  = self.bf
        tf = self.tf
        tw = self.tw
        d  = self.d

        # Shear from Saritas and Filippou (2009)
        # Ratio of total flange area to web area
        alpha = 2*b*tf/d/(2*tw)
        # NOTE: This is 1/beta_S where beta_S is Afsin's beta
        beta = (1+3*alpha)*(2/3)/((1+2*alpha)**2-2/3*(1+2*alpha)+1/5)
        def psi(y, z):
            # webs
            if abs(y) < (d/2-tf):
                return 0 #beta*((1+2*alpha) - (2*y/d)**2) - 1 #+ 1
            # flange
            else:
                return 0 #beta*(2*alpha)*(z/b) - 1

        mesh = self._create_mesh(mesh_scale=mesh_scale)

        return TriangleSection.from_meshio(mesh, warp_shear=psi)



class Rectangle(_BasicShape):
    def __init__(self, b, d, **kwds):
        self.b = b
        self.d = d
        # using _r_ prefix to avoid clobbering parent variables,
        # but this may not actually be necessary
        super().__init__(mesh_size=min(b,d)/10, **kwds)

    def exterior(self):
        b = self.b
        d = self.d
        return np.array([
            [-b / 2,  -d / 2],
            [ b / 2,  -d / 2],
            [ b / 2,   d / 2],
            [-b / 2,   d / 2],
        ])


class HollowRectangle(_BasicShape):
    def __init__(self, b, d, tf, tw, t=None, **kwds):
        self.b = b
        self.d = d
        self.t = t

        super().__init__(mesh_size=min(tf, tw)/10, **kwds)

    def exterior(self):
        b = self.b
        d = self.d
        return np.array([
            [-b / 2, -d / 2],
            [ b / 2, -d / 2],
            [ b / 2,  d / 2],
            [-b / 2,  d / 2],
        ], dtype=float)
    
    def interior(self):
        b = self.b
        d = self.d
        tf, tw = self.tf, self.tw
        return [np.array([
            [-(b/2 - tw), -(d/2 - tf)],
            [ (b/2 - tw), -(d/2 - tf)],
            [ (b/2 - tw),  (d/2 - tf)],
            [-(b/2 - tw),  (d/2 - tf)],
        ], dtype=float)]

    def _create_patches(self, mesh_scale=None):
        if mesh_scale is None:
            mesh_scale = 1/5
        t = self.t
        x1 = self.b/2 - t
        x2 = self.b/2
        y1 = self.d/2 - t
        y2 = self.d/2


        return create_mesh(mesh_size=t*mesh_scale, patches=[
            patch.rect(corners=[[-x2, -y2], [ x2, -y1]]),
            patch.rect(corners=[[-x2, -y1], [-x1,  y1]]),
            patch.rect(corners=[[ x1, -y1], [ x2,  y1]]),
            patch.rect(corners=[[-x2,  y1], [ x2,  y2]]),
        ])


class Channel(_BasicShape):
    """
    _  ___________
      |__________|
      | |
      | |
      |o|
      | |
      | |_________
      |__________|



    """
    def __init__(self, d, b, tf, tw=None, **kwds):
        self.tf = tf
        self.tw = tw if tw is not None else tf
        self.d = d
        self.b = b

        super().__init__(mesh_size=min(tf, tw)/4, **kwds)

    def exterior(self):
        d = self.d
        b = self.b
        t = self.tf
        w = self.tw

        y_top = d / 2.0
        y_bot = -d / 2.0

        pts = np.array([
            [-w / 2, y_top],            # 1  top-left
            [ b - w / 2, y_top],        # 2  top-right
            [ b - w / 2, y_top - t],    # 3  down into flange
            [  w / 2, y_top - t],       # 4  over to web
            [  w / 2, y_bot + t],       # 5  down web
            [ b - w / 2, y_bot + t],    # 6  out to bottom flange
            [ b - w / 2, y_bot],        # 7  bottom-right
            [-w / 2, y_bot],            # 8  bottom-left
        ], dtype=float)

        return pts



class Angle(_BasicShape):
    def __init__(self, t, b, d, **kwds):
        self.t = t
        self.b = b
        self.d = d
        super().__init__(mesh_size=min(t, b)/4, **kwds)

    def exterior(self):
        t = self.t
        b = self.b
        d = self.d

        pts = np.array([
            [ b - t / 2,  t / 2],          # 1  outer top-right
            [-t / 2,      t / 2],          # 2  outer top-left
            [-t / 2, -d + t / 2],          # 3  outer bottom-left
            [ t / 2, -d + t / 2],          # 4  inner bottom-left of vertical leg
            [ t / 2,     -t / 2],          # 5  inner corner
            [ b - t / 2, -t / 2],          # 6  outer bottom-right of horizontal leg
        ], dtype=float)

        return pts

