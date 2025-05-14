from ..utility import find_row, find_rows, UnimplementedInstance
import numpy as np
import warnings
from veux.frame import SectionGeometry

_CIRCLE_DIVS = 40

def _HatSO3(vec):
    """Construct a skew-symmetric matrix from a 3-vector."""
    return np.array([
        [0, -vec[2], vec[1]],
        [vec[2], 0, -vec[0]],
        [-vec[1], vec[0], 0]
    ])

def _ExpSO3(vec):
    """
    Exponential map for SO(3).
    Satisfies ExpSO3(vec) == expm(skew(vec)).
    """
    vec = np.asarray(vec)
    if vec.shape != (3,):
        raise ValueError("Input must be a 3-vector.")

    theta = np.linalg.norm(vec)
    if theta < 1e-8:  # Small-angle approximation
        return np.eye(3) + _HatSO3(vec) + 0.5 * (_HatSO3(vec) @ _HatSO3(vec))
    else:
        K = _HatSO3(vec / theta)  # Normalized skew matrix
        return np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)



def create_frame_sections(csi, model, conv):
    for sect in csi.get("FRAME SECTION PROPERTIES 01 - GENERAL", []):

        if not conv.identify("AnalSect", "section",     sect["SectionName"]) and \
           not conv.identify("AnalSect", "integration", sect["SectionName"]):

            if (s:= _create_section(csi, sect, model, conv)) is not None:
                continue

            if (s := _create_integration(csi, sect, model, conv)) is not None:
                continue

            assert False, sect
            conv.log(UnimplementedInstance(f"FrameSection.Shape={sect['Shape']}"))
    return

def collect_geometry(csi, elem_maps=None, conv=None):
    """
    collect section geometry
    """

    frame_types = {
        row["SectionName"]: FrameQuadrature.from_table(csi, row)
        for row in csi.get("FRAME SECTION PROPERTIES 01 - GENERAL", [])
    }

    frame_assigns = {}
    for row in csi.get("FRAME SECTION ASSIGNMENTS",[]):

        if row["MatProp"] != "Default":
            if conv is not None:
                conv.log(UnimplementedInstance("FrameSection.MatProp", row["MatProp"]))
            else:
                warnings.warn(f"Material property {row['MatProp']} not implemented.")

        if row["AnalSect"] in frame_types and frame_types[row["AnalSect"]] is not None:
            if frame_types[row["AnalSect"]].geometry() is None:
                warnings.warn(f"No geometry for {row['AnalSect']}")
                continue
            frame_assigns[row["Frame"]] = frame_types[row["AnalSect"]].geometry()
            


    # Skew angles
    E2 = np.array([0, 0,  1])
    for frame in frame_assigns:
        skew_assign = find_row(csi.get("FRAME END SKEW ANGLE ASSIGNMENTS", []),
                        Frame=frame)
        
        if skew_assign: #and skew["SkewI"] != 0 and skew["SkewJ"] != 0: # and len(frame_assigns[frame].shape) == 2
            for i,skew in zip((0,-1), ("SkewI", "SkewJ")):
                exterior = frame_assigns[frame][i].exterior()
                interior = frame_assigns[frame][i].interior()

                R = _ExpSO3(skew_assign[skew]*np.pi/180*E2)
                frame_assigns[frame][i] = SectionGeometry(interior=[np.array([[(R@point)[0], *point[1:]] for point in hole]) for hole in interior],
                                                          exterior=np.array([[(R@point)[0], *point[1:]]  for point in exterior])
                )

    if elem_maps is not None:
        return {
            elem_maps.get(name,name): val for name, val in frame_assigns.items()
        }
    else:
        return frame_assigns
    

def section_geometry(csi, prop_01):
    if isinstance(prop_01, str):
        name = prop_01
        prop_01 = find_row(csi.get("FRAME SECTION PROPERTIES 01 - GENERAL",[]), SectionName=name)
        if prop_01 is None:
            prop_01 = find_row(csi.get("FRAME SECTION PROPERTIES - BRIDGE OBJECT FLAGS",[]), SectionName=name)
            if prop_01 is None:
                raise ValueError(f"Section {name} not found in either table.")
    else:
        name = prop_01["SectionName"]

    exterior = None
    interior = []
    if "Shape" not in prop_01:
        return 

    if prop_01["Shape"] == "Circle":
        r = prop_01["t3"]/2
        exterior = np.array([
            [np.sin(x)*r, np.cos(x)*r] for x in np.linspace(0, np.pi*2, _CIRCLE_DIVS)
        ])
    elif prop_01["Shape"] == "Rectangular":
        # TODO: Check if 2/3 axes are correct
        exterior = np.array([
            [ prop_01["t2"]/2,   prop_01["t3"]/2],
            [ prop_01["t2"]/2,  -prop_01["t3"]/2],
            [-prop_01["t2"]/2,  -prop_01["t3"]/2],
            [-prop_01["t2"]/2,   prop_01["t3"]/2],
            [ prop_01["t2"]/2,   prop_01["t3"]/2],
        ]) + np.array([prop_01.get("CGOffset2",0), prop_01.get("CGOffset3", 0)])

    elif prop_01["Shape"] == "SD Section":
        prop_sd = find_row(csi.get("SECTION DESIGNER PROPERTIES 01 - GENERAL", []), SectionName=name)
        if prop_sd["nCaltransCr"] == 1:
            circle = find_row(csi.get("SECTION DESIGNER PROPERTIES 24 - SHAPE CALTRANS CIRCLE", []), SectionName=name)
            assert circle["Height"] == circle["Width"]
            r = circle["Height"]/2
            exterior = np.array([
                [np.sin(x)*r, np.cos(x)*r] for x in np.linspace(0, np.pi*2, _CIRCLE_DIVS)
            ])

        elif prop_sd["nPolygon"] > 0:
            if prop_sd["nPolygon"] != prop_sd["nTotalShp"]:
                # unimplemented
                return
            polygon_data = csi.get("SECTION DESIGNER PROPERTIES 16 - SHAPE POLYGON", [])

            exterior =  np.array([
                [row["X"], row["Y"]]
                for row in find_rows(polygon_data, SectionName = name) if row.get("ShapeName","")=="Polygon1"
            ])
            if len(exterior) == 0:
                return
            
            for hole in find_rows(polygon_data, SectionName = name, ShapeMat="Opening"):
                interior.append(np.array([
                    [row["X"], row["Y"]]
                    for row in find_rows(polygon_data, ShapeName=hole["ShapeName"])
                ]))
        else:
            warnings.warn(f"Unimplemented section designer section.")

    elif prop_01["Shape"] == "Bridge Section":
        polygon_data = csi.get("FRAME SECTION PROPERTIES 06 - POLYGON DATA", [])
        exterior_row = find_row(polygon_data, SectionName = name, Opening=False)
        exterior =  np.array([
            [row["X"], row["Y"]]
            for row in find_rows(polygon_data, SectionName = name) if row["Polygon"] == exterior_row["Polygon"]
        ])
        ref = (exterior_row["RefPtX"], exterior_row["RefPtY"])

        for i in range(len(exterior)):
            exterior[i] -= ref


        for hole in find_rows(polygon_data, SectionName = name, Opening=True):
            interior.append(np.array([
                [row["X"], row["Y"]]
                for row in find_rows(polygon_data, Polygon=hole["Polygon"])
            ]))
            for i in range(len(interior[-1])):
                interior[-1][i] -= ref


    if exterior is not None:
        return SectionGeometry(exterior, interior=interior)




def section_mesh(csi, prop_01, engine=None):

    from shps.frame.mesh import sect2meshpy
    shape = section_geometry(csi, prop_01)
    if engine is None:
        shape = (
            shape.exterior(plane=True),
            shape.interior(plane=True)
        )
        return sect2meshpy(shape, 0.5)

    geometry = section_geometry(csi, prop_01)
    exterior = geometry.exterior(plane=True)
    interior = geometry.interior(plane=True)

    import gmsh
    import meshio
    gmsh.initialize()
    gmsh.model.add("section")
    # Add exterior points
    exterior_points = [gmsh.model.geo.addPoint(x, y, 0) for x, y in exterior]
    exterior_loop = gmsh.model.geo.addCurveLoop([gmsh.model.geo.addSpline(exterior_points)])
    
    # Add interior points (holes)
    interior_loops = []
    for hole in interior:
        hole_points = [gmsh.model.geo.addPoint(x, y, 0) for x, y in hole]
        interior_loops.append(gmsh.model.geo.addCurveLoop([gmsh.model.geo.addSpline(hole_points)]))
    
    # Create plane surface with holes
    surface = gmsh.model.geo.addPlaneSurface([exterior_loop] + interior_loops)
    
    # Synchronize to create the surface
    gmsh.model.geo.synchronize()
    
    # Generate 2D mesh
    gmsh.model.mesh.generate(2)
    
    # Extract mesh data
    nodes = gmsh.model.mesh.getNodes()
    elements = gmsh.model.mesh.getElements()
    
    # Convert to meshio format
    mesh = meshio.Mesh(
        points=nodes[1].reshape(-1, 3)[:, :2],  # Only take x, y coordinates
        cells={"triangle": elements[2][0].reshape(-1, 3) - 1}  # Convert to 0-based indexing
    )
    
    # Finalize gmsh
    gmsh.finalize()
    
    return mesh
    exterior_loop = gmsh.model.geo.addCurveLoop([
        gmsh.model.geo.addSpline([gmsh.model.geo.addPoint(x, y, 0) for x, y in exterior])
    ])
    
    # Add interior points (holes)
    interior_loops = []
    for hole in interior:
        interior_loops.append(gmsh.model.geo.addCurveLoop([gmsh.model.geo.addSpline([gmsh.model.geo.addPoint(x, y, 0) for x, y in hole])]))
    
    # Create plane surface with holes
    surface = gmsh.model.geo.addPlaneSurface([exterior_loop] + interior_loops)
    
    # Synchronize to create the surface
    gmsh.model.geo.synchronize()
    
    # Generate 2D mesh
    gmsh.model.mesh.generate(2)
    
    # Extract mesh data
    nodes = gmsh.model.mesh.getNodes()
    elements = gmsh.model.mesh.getElements()
    
    # Convert to meshio format
    mesh = meshio.Mesh(
        points=nodes[1].reshape(-1, 3)[:, :2],  # Only take x, y coordinates
        cells={"triangle": elements[2][0].reshape(-1, 3) - 1}  # Convert to 0-based indexing
    )
    
    # Finalize gmsh
    gmsh.finalize()
    
    return mesh



class _FrameSection:
    def __init__(self, geometry=None):
        self._geometry = geometry

    @classmethod
    def from_table(cls, prop_01, csi):
        s = _FrameSection()
        s._prop_01 = prop_01
        s._csi = csi
        return s

    # def geometry(self):
    #     return [i for i in self._geometry]

    def add_to(self, model, conv, name=None):
        prop_01 = self._prop_01
        csi = self._csi

        material = find_row(csi.get("MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES", []),
                            Material=prop_01["Material"]
        )

        if material is None:
            print(prop_01)

        if "G12" in material:
            model.section("FrameElastic",
                            conv.define("AnalSect", "section", name), #self.index,
                            A  = prop_01["Area"],
                            Ay = prop_01["AS3"],
                            Az = prop_01["AS2"],
                            Iz = prop_01["I33"],
                            Iy = prop_01["I22"],
                            J  = prop_01["TorsConst"],
                            E  = material["E1"],
                            G  = material["G12"]
            )

    def cnn(self):
        pass 

    def cmm(self):
        pass


class FrameQuadrature:
    def __init__(self, sections, locations=None, geometry: list=None):
        self._sections  = sections
        self._locations = locations
        self._geometry  = geometry

    @classmethod
    def from_table(cls, csi, prop_01):
        # 1)
        if prop_01["Shape"] != "Nonprismatic":
            geometry = section_geometry(csi, prop_01)
            if geometry is not None:
                geometry = [geometry, geometry]
            section = _FrameSection.from_table(csi, prop_01)
            return FrameQuadrature([section, section],
                                   geometry=geometry)


        row = find_row(csi.get("FRAME SECTION PROPERTIES 05 - NONPRISMATIC", []),
                        SectionName=prop_01["SectionName"])

        # 2)
        if row["StartSect"] == row["EndSect"]:
            si = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"], 
                                SectionName=row["StartSect"])

            assert si is not None

            geometry = section_geometry(csi, si)
            if geometry is not None:
                geometry = [geometry, geometry]
            section = _FrameSection.from_table(csi, prop_01)
            return FrameQuadrature([section, section], geometry=geometry)

        # 3)
        else:
            si = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
                            SectionName=row["StartSect"])
            sj   = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
                            SectionName=row["EndSect"])

            if si["Shape"] == sj["Shape"] and si["Shape"] in {"Circle"}:
                circumference = np.linspace(0, np.pi*2, _CIRCLE_DIVS)
                exteriors = np.array([
                    [[np.sin(x)*r, np.cos(x)*r] for x in circumference]
                    for r in np.linspace(si["t3"]/2, sj["t3"]/2, 2)
                ])
                sections = []
                return FrameQuadrature(sections, 
                                       geometry = [SectionGeometry(exterior) for exterior in exteriors])

    def sections(self):
        return self._sections

    def locations(self):
        pass 

    def geometry(self):
        if self._geometry is not None:
            return [i for i in self._geometry]

    def weights(self):
        pass


def _create_section(csi, prop_01, model, conv):

    # 1)
    if prop_01["Shape"] not in {"Nonprismatic"}:
        s = _FrameSection.from_table(prop_01, csi)
        s.add_to(model, conv, name=prop_01["SectionName"])
        return s


    segments = find_rows(csi.get("FRAME SECTION PROPERTIES 05 - NONPRISMATIC",[]),
                            SectionName=prop_01["SectionName"])

    # 2)
    if prop_01["Shape"] == "Nonprismatic" and \
            len(segments) != 1: #section["NPSecType"] == "Advanced":

        # TODO: Currently just treating advanced as normal prismatic section

        assert all(segment["StartSect"] == segment["EndSect"] for segment in segments)

        if not conv.identify("AnalSect", "section", prop_01["SectionName"]) : #segments[0]["StartSect"]):
            # find properties
            p = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
                                SectionName=segments[0]["StartSect"]
            )
            assert p is not None
            s = _FrameSection.from_table(p, csi)
            s.add_to(model, conv, name=prop_01["SectionName"]) #segments[0]["StartSect"])
            return s



def _create_integration(csi, prop_01, model, conv):
    # 3)
    segments = find_rows(csi["FRAME SECTION PROPERTIES 05 - NONPRISMATIC"],
                            SectionName=prop_01["SectionName"])
    if prop_01["Shape"] != "Nonprismatic" or len(segments) != 1: 
        return None
    #assert section["NPSecType"] == "Default":
    assert len(segments) == 1

    segment = segments[0]

    # Create property interpolation
    def interpolate(point, prop):
        si = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
                            SectionName=segment["StartSect"]
        )
        sj = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
                            SectionName=segment["EndSect"]
        )
        # TODO: Taking material from first section assumes si and sj have the same
        # material
        material = find_row(csi["MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES"],
                            Material=si["Material"]
        )

        if prop in material:
            start = end = material[prop]
        else:
            start = si[prop]
            end = sj[prop]

        power = {
                "Linear":    1,
                "Parabolic": 2,
                "Cubic":     3
        }[segment.get(f"E{prop}Var", "Linear")]

        return start*(1 + point*((end/start)**(1/power)-1))**power


    # Define a numerical integration scheme

    from numpy.polynomial.legendre import leggauss
    nip = 5
    sections = []
    for x,wi in zip(*leggauss(nip)):
        xi = (1+x)/2
        #tag = self.index+off
        tag = conv.define("AnalSect", "section")

        model.section("FrameElastic", 
                        tag, #self.index+off, #, #
                        A  = interpolate(xi, "Area"),
                        Ay = interpolate(xi, "AS2"),
                        Az = interpolate(xi, "AS2"),
                        Iz = interpolate(xi, "I33"),
                        Iy = interpolate(xi, "I22"),
                        J  = interpolate(xi, "TorsConst"),
                        E  = interpolate(xi, "E1"),
                        G  = interpolate(xi, "G12")
        )


        # self.integration.append((tag, xi, wi/2))
        sections.append((tag, xi, wi/2))

    model.beamIntegration("UserDefined",
                            conv.define("AnalSect", "integration", prop_01["SectionName"]),
                            len(sections),
                            tuple(i[0] for i in sections),
                            tuple(i[1] for i in sections),
                            tuple(i[2] for i in sections))

    return sections
