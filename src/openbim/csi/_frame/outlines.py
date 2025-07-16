#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
from ..utility import find_row, find_rows
import numpy as np
import warnings
from veux.frame import SectionGeometry

from .section import CIRCLE_DIVS as _CIRCLE_DIVS, _CsiSection


def _find_polygons(csi, prop_01) -> SectionGeometry:
    """
    was previously section_geometry(csi, prop_01)

    POST-PROCESSING
    """
    #
    if isinstance(prop_01, str):
        name = prop_01
        prop_01 = find_row(csi.get("FRAME SECTION PROPERTIES 01 - GENERAL",[]), SectionName=name)
        if prop_01 is None:
            prop_01 = find_row(csi.get("FRAME SECTION PROPERTIES - BRIDGE OBJECT FLAGS",[]), SectionName=name)
            if prop_01 is None:
                raise ValueError(f"Section {name} not found in either table.")
    else:
        name = prop_01["SectionName"]

    s = _CsiSection(csi, prop_01)._create_model(render_only=True)
    if s is not None:
        return SectionGeometry(s.exterior(), interior=s.interior())


# def _section_mesh(csi, prop_01, engine=None):
#     """
#     """
#     from shps.frame.mesh import sect2meshpy
#     shape = _find_polygons(csi, prop_01)
#     if engine is None:
#         shape = (
#             shape.exterior(plane=True),
#             shape.interior(plane=True)
#         )
#         return sect2meshpy(shape, 0.5)

#     geometry = _find_polygons(csi, prop_01)
#     exterior = geometry.exterior(plane=True)
#     interior = geometry.interior(plane=True)

#     import gmsh
#     import meshio
#     gmsh.initialize()
#     gmsh.model.add("section")
#     # Add exterior points
#     exterior_points = [gmsh.model.geo.addPoint(x, y, 0) for x, y in exterior]
#     exterior_loop = gmsh.model.geo.addCurveLoop([gmsh.model.geo.addSpline(exterior_points)])

#     # Add interior points (holes)
#     interior_loops = []
#     for hole in interior:
#         hole_points = [gmsh.model.geo.addPoint(x, y, 0) for x, y in hole]
#         interior_loops.append(gmsh.model.geo.addCurveLoop([gmsh.model.geo.addSpline(hole_points)]))

#     # Create plane surface with holes
#     surface = gmsh.model.geo.addPlaneSurface([exterior_loop] + interior_loops)

#     # Synchronize to create the surface
#     gmsh.model.geo.synchronize()

#     # Generate 2D mesh
#     gmsh.model.mesh.generate(2)
    
#     # Extract mesh data
#     nodes = gmsh.model.mesh.getNodes()
#     elements = gmsh.model.mesh.getElements()
    
#     # Convert to meshio format
#     mesh = meshio.Mesh(
#         points=nodes[1].reshape(-1, 3)[:, :2],  # Only take x, y coordinates
#         cells={"triangle": elements[2][0].reshape(-1, 3) - 1}  # Convert to 0-based indexing
#     )
    
#     # Finalize gmsh
#     gmsh.finalize()
    
#     return mesh


class FrameQuadrature:
    def __init__(self, sections, locations=None, geometry: list=None, name=None):
        self._sections  = sections
        self._locations = locations
        self._geometry  = geometry
        self._name = name 

    def __repr__(self):
        if self._name is not None:
            return f"FrameQuadrature({self._name})"
        else:
            return f"FrameQuadrature({len(self._sections)} sections, {len(self._geometry)} geometries)"

    @classmethod
    def from_table(cls, csi, prop_01):
        # 1)
        name = prop_01["SectionName"]
        if prop_01["Shape"] != "Nonprismatic":
            geometry = _find_polygons(csi, prop_01)
            if geometry is not None:
                geometry = [geometry, geometry]

            # section = _FrameSection.from_table(csi, prop_01)

            return FrameQuadrature([], #[section, section],
                                   geometry=geometry, name=name)

        row = find_row(csi.get("FRAME SECTION PROPERTIES 05 - NONPRISMATIC", []),
                        SectionName=prop_01["SectionName"])

        # 2)
        if row["StartSect"] == row["EndSect"]:
            si = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"], 
                                SectionName=row["StartSect"])

            assert si is not None

            geometry = _find_polygons(csi, si)
            if geometry is not None:
                geometry = [geometry, geometry]
            # section = _FrameSection.from_table(csi, prop_01)
            return FrameQuadrature([], #[section, section],
                                   geometry=geometry, name=name)

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
                return FrameQuadrature([], 
                                       geometry = [SectionGeometry(exterior) for exterior in exteriors])

    # def sections(self):
    #     return self._sections

    def locations(self):
        pass 

    def geometry(self):
        if self._geometry is not None:
            return [i for i in self._geometry]

    def weights(self):
        pass


def collect_geometry(csi, elem_maps=None):
    """
    collect a mapping of SectionGeometry for all frames elements in a model.

    This is a POST-PROCESSING function, used for rendering extruded frames.

    TODO: Reimplement as loop over iter_sections()
    """

    frame_types = {
        row["SectionName"]: FrameQuadrature.from_table(csi, row)
        for row in csi.get("FRAME SECTION PROPERTIES 01 - GENERAL", [])
    }

    frame_assigns = {}
    for row in csi.get("FRAME SECTION ASSIGNMENTS",[]):
        if row["MatProp"] != "Default":
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

