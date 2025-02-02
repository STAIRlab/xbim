from .utility import find_row, find_rows, UnimplementedInstance
import numpy as np
import warnings
from veux.frame import SectionGeometry
from collections import defaultdict

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


def collect_sections(csi, elem_maps=None, conv=None):

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
            frame_assigns[row["Frame"]] = [s.geometry() for s in frame_types[row["AnalSect"]].sections()]
            # np.array([
            #     section.geometry().exterior() for section in
            #     frame_types[row["AnalSect"]].sections()
            # ])


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
                                                          exterior=np.array([[(R@point)[0], *point[1:]] for point in exterior])
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
        prop_01 = find_row(csi.get("FRAME SECTION PROPERTIES 01 - GENERAL"), SectionName=name)
    else:
        name = prop_01["SectionName"]

    exterior = None
    interior = []

    if "Shape" not in prop_01:
        return 

    if prop_01["Shape"] == "Circle":
        r = prop_01["t3"]/2
        exterior = np.array([
            [np.sin(x)*r, np.cos(x)*r] for x in np.linspace(0, np.pi*2, 40)
        ])

    elif prop_01["Shape"] == "SD Section":
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


class _FrameSection:
    def __init__(self, geometry):
        self._geometry = geometry

    def geometry(self):
        return self._geometry

    def add_to(self, model, type, tag):
        pass

    def cnn(self):
        pass 

    def cmm(self):
        pass


class FrameQuadrature:
    def __init__(self, sections):
        self._sections = sections

    @classmethod
    def from_table(cls, csi, prop_01):

        if prop_01["Shape"] != "Nonprismatic":
            geometry = section_geometry(csi, prop_01)
            if geometry is None:
                return None
            section = _FrameSection(
                geometry=geometry
            )
            return FrameQuadrature([section, section])

        else:
            row = find_row(csi.get("FRAME SECTION PROPERTIES 05 - NONPRISMATIC", []),
                           SectionName=prop_01["SectionName"])

            if row["StartSect"] == row["EndSect"]:
                start = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"], 
                                 SectionName=row["StartSect"])

                if start is None:
                    print(f"Section {row['StartSect']} not found.")

                section = _FrameSection(
                    geometry=section_geometry(csi, start)
                )
                return FrameQuadrature([section, section])

            start = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
                            SectionName=row["StartSect"])
            end   = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
                            SectionName=row["EndSect"])

            if start["Shape"] == end["Shape"] and start["Shape"] in {"Circle"}:
                circumference = np.linspace(0, np.pi*2, 40)
                exteriors = np.array([
                    [[np.sin(x)*r, np.cos(x)*r] for x in circumference]
                    for r in np.linspace(start["t3"]/2, end["t3"]/2, 2)
                ])

                return FrameQuadrature([
                    _FrameSection(geometry=SectionGeometry(exterior)) for exterior in exteriors
                ])


    def sections(self):
        return self._sections

    def locations(self):
        pass 

    def weights(self):
        pass

