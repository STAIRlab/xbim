import warnings 
from .utility import find_row, find_rows
import numpy as np

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


class _Section:
    def __init__(self, name: str, csi: dict,
                 index: int, model, library, conv):
        self.index = index
        self.name = name
        self.integration = []

        self._create(csi, model, library, conv)

    def _create(self, csi, model, library, conv):
        pass


class _ShellSection(_Section):
    def _create(self, csi, model, library, conv):

        section = find_row(csi["AREA SECTION PROPERTIES"],
                           Section=self.name
        )

        if section is None:
            print(self.name)

        material = find_row(csi["MATERIAL PROPERTIES 01 - GENERAL"],
                            Material=section["Material"]
        )

        material = find_row(csi["MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES"],
                            Material=section["Material"]
        )
        model.section("ElasticMembranePlateSection", self.index,
                      material["E1"],  # E
                      material["G12"]/(2*material["E1"]) - 1, # nu
                      section["Thickness"],
                      material["UnitMass"]
        )
        self.integration.append(self.index)


# def collect_outlines(csi, elem_maps=None):
#     polygon_data = csi.get("FRAME SECTION PROPERTIES 06 - POLYGON DATA", 
#                            csi.get("SECTION DESIGNER PROPERTIES 16 - SHAPE POLYGON", 
#                            []))

#     names = set(
#             row["SectionName"] for row in polygon_data
#     )

#     exteriors = {
#             s: np.array([
#                 [row["X"], row["Y"]]
#                 for row in polygon_data 
#                 if (row.get("Polygon", 0) == 1 or row.get("ShapeName","")=="Polygon1") and row["SectionName"] == s
#             ])
#             for s in names
#     }

#     # Adjust from reference point
#     for row in polygon_data:
#         if row.get("Polygon", 0) == 1 and row["Point"] == 1:
#             ref = (row["RefPtX"], row["RefPtY"])

#             for i in range(len(exteriors[row["SectionName"]])):
#                 exteriors[row["SectionName"]][i] -= ref


#     for row in csi.get("FRAME SECTION PROPERTIES 01 - GENERAL", []):
#         if "Shape" in row and row["Shape"] == "Circle":
#             r = row["t3"]/2
#             exteriors[row["SectionName"]] = np.array([
#                 [np.sin(x)*r, np.cos(x)*r] for x in np.linspace(0, np.pi*2, 40)
#             ])

#     for row in csi.get("FRAME SECTION PROPERTIES 05 - NONPRISMATIC", []):
#         if row["StartSect"] == row["EndSect"] and row["StartSect"] in exteriors:
#             exteriors[row["SectionName"]] = exteriors[row["StartSect"]]
#         else:
#             start = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
#                              SectionName=row["StartSect"])
#             end   = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
#                              SectionName=row["EndSect"])

#             if start["Shape"] == end["Shape"] and start["Shape"] in {
#                     "Circle"
#                 }:
#                 exteriors[row["SectionName"]] = np.array([
#                     [[0, np.sin(x)*r, np.cos(x)*r] for x in np.linspace(0, np.pi*2, 40)]
#                     for r in np.linspace(start["t3"]/2, end["t3"]/2, 2)
#                 ])


#     frame_sections = {
#             row["Frame"]: exteriors[row["AnalSect"]]
#             for row in csi.get("FRAME SECTION ASSIGNMENTS",[]) if row["AnalSect"] in exteriors
#     }

#     # Skew angles
#     for frame in frame_sections:
#         skew = find_row(csi.get("FRAME END SKEW ANGLE ASSIGNMENTS", []),
#                         Frame=frame)

#         if skew and len(frame_sections[frame].shape) == 2: #and skew["SkewI"] != 0 and skew["SkewJ"] != 0:
#             E2 = np.array([0, 0,  1])
#             RI = ExpSO3(skew["SkewI"]*np.pi/180*E2)
#             RJ = ExpSO3(skew["SkewJ"]*np.pi/180*E2)
#             frame_sections[frame] = np.array([
#                 [[(RI@[0, *point])[0], *point] for point in frame_sections[frame]],
#                 [[(RJ@[0, *point])[0], *point] for point in frame_sections[frame]],
#             ])

#     if elem_maps is not None:
#         return {
#             elem_maps.get(name,name): val for name, val in frame_sections.items()
#         }
#     else:
#         return frame_sections


def create_shell_sections(csi, model, library, conv):
    for assign in csi.get("AREA SECTION ASSIGNMENTS", []):
        if assign["Section"] not in library["shell_sections"]:
            library["shell_sections"][assign["Section"]] = \
              _ShellSection(assign["Section"], csi, tag, model, library)
            tag += len(library["shell_sections"][assign["Section"]].integration)
