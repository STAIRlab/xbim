import warnings 
from .utility import find_row, find_rows
import numpy as np

def skew(vec):
    """Construct a skew-symmetric matrix from a 3-vector."""
    return np.array([
        [0, -vec[2], vec[1]],
        [vec[2], 0, -vec[0]],
        [-vec[1], vec[0], 0]
    ])

def ExpSO3(vec):
    """
    Exponential map for SO(3).
    Satisfies ExpSO3(vec) == expm(skew(vec)).
    """
    vec = np.asarray(vec)
    if vec.shape != (3,):
        raise ValueError("Input must be a 3-vector.")

    theta = np.linalg.norm(vec)
    if theta < 1e-8:  # Small-angle approximation
        return np.eye(3) + skew(vec) + 0.5 * (skew(vec) @ skew(vec))
    else:
        K = skew(vec / theta)  # Normalized skew matrix
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


def collect_outlines(csi, elem_maps=None):
    polygon_data = csi.get("FRAME SECTION PROPERTIES 06 - POLYGON DATA", 
                           csi.get("SECTION DESIGNER PROPERTIES 16 - SHAPE POLYGON", 
                           []))

    names = set(
            row["SectionName"] for row in polygon_data
    )

    polygons = {
            s: np.array([
                [row["X"], row["Y"]]
                for row in polygon_data if (row.get("Polygon", 0) == 1 or row.get("ShapeName","")=="Polygon1") and row["SectionName"] == s
            ])
            for s in names
    }

    # Adjust from reference point
    for row in polygon_data:
        if row.get("Polygon", 0) == 1 and row["Point"] == 1:
            ref = (row["RefPtX"], row["RefPtY"])

            for i in range(len(polygons[row["SectionName"]])):
                polygons[row["SectionName"]][i] -= ref


    for row in csi.get("FRAME SECTION PROPERTIES 01 - GENERAL", []):
        if "Shape" in row and row["Shape"] == "Circle":
            r = row["t3"]/2
            polygons[row["SectionName"]] = np.array([
                [np.sin(x)*r, np.cos(x)*r] for x in np.linspace(0, np.pi*2, 40)
            ])

    for row in csi.get("FRAME SECTION PROPERTIES 05 - NONPRISMATIC", []):
        if row["StartSect"] == row["EndSect"] and row["StartSect"] in polygons:
            polygons[row["SectionName"]] = polygons[row["StartSect"]]
        else:
            start = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
                             SectionName=row["StartSect"])
            end   = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
                             SectionName=row["EndSect"])

            if start["Shape"] == end["Shape"] and start["Shape"] in {
                    "Circle"
                }:
                polygons[row["SectionName"]] = np.array([
                    [[0, np.sin(x)*r, np.cos(x)*r] for x in np.linspace(0, np.pi*2, 40)]
                    for r in np.linspace(start["t3"]/2, end["t3"]/2, 2)
                ])


    frame_sections = {
            row["Frame"]: polygons[row["AnalSect"]]
            for row in csi.get("FRAME SECTION ASSIGNMENTS",[]) if row["AnalSect"] in polygons
    }

    # Skew angles
    for frame in frame_sections:
        skew = find_row(csi.get("FRAME END SKEW ANGLE ASSIGNMENTS", []),
                        Frame=frame)

        if skew and len(frame_sections[frame].shape) == 2: #and skew["SkewI"] != 0 and skew["SkewJ"] != 0:
            E2 = np.array([0, 0,  1])
            RI = ExpSO3(skew["SkewI"]*np.pi/180*E2)
            RJ = ExpSO3(skew["SkewJ"]*np.pi/180*E2)
            frame_sections[frame] = np.array([
                [[(RI@[0, *point])[0], *point] for point in frame_sections[frame]],
                [[(RJ@[0, *point])[0], *point] for point in frame_sections[frame]],
            ])

    if elem_maps is not None:
        return {
            elem_maps.get(name,name): val for name, val in frame_sections.items()
        }
    else:
        return frame_sections



def create_frame_sections(csi, model, library, conv):
    tag = 1
    for assign in csi.get("FRAME SECTION ASSIGNMENTS", []):

        if assign["AnalSect"] not in library["frame_sections"]:

            library["frame_sections"][assign["AnalSect"]] = \
              _FrameSection(assign["AnalSect"], csi, tag, model, library, conv)

            tag += len(library["frame_sections"][assign["AnalSect"]].integration)

    return tag


def create_shell_sections(csi, model, library, conv):
    for assign in csi.get("AREA SECTION ASSIGNMENTS", []):
        if assign["Section"] not in library["shell_sections"]:
            library["shell_sections"][assign["Section"]] = \
              _ShellSection(assign["Section"], csi, tag, model, library)
            tag += len(library["shell_sections"][assign["Section"]].integration)


class FrameSection:
    def __init__(self, name, csi):
        pass

    def exterior(self, name, csi):
        prop_01 = find_row(csi.get("FRAME SECTION PROPERTIES 01 - GENERAL", []),
                           SectionName=name)

        if "Shape" not in prop_01:
            return 
        
        if prop_01["Shape"] == "Circle":
            r = prop_01["t3"]/2
            exterior = np.array([
                [np.sin(x)*r, np.cos(x)*r] for x in np.linspace(0, np.pi*2, 40)
            ])

        elif prop_01["Shape"] in {"SD Section", "Bridge Section"}:
            polygon_data = csi.get("FRAME SECTION PROPERTIES 06 - POLYGON DATA", 
                            csi.get("SECTION DESIGNER PROPERTIES 16 - SHAPE POLYGON", 
                            []))
            
            exterior =  np.array([
                [row["X"], row["Y"]]
                for row in polygon_data if row["SectionName"] == name and (row.get("Polygon", 0) == 1 or row.get("ShapeName","")=="Polygon1")
            ])
            

        elif prop_01["Shape"] == "Nonprismatic":
            row = find_row(csi.get("FRAME SECTION PROPERTIES 05 - NONPRISMATIC", []),
                           SectionName=name)

            if row["StartSect"] == row["EndSect"]:
                exterior = FrameSection().exterior(row["StartSect"], csi)
                if exterior is not None:
                    return exterior

            start = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
                            SectionName=row["StartSect"])
            end   = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
                            SectionName=row["EndSect"])

            if start["Shape"] == end["Shape"] and start["Shape"] in {
                    "Circle"
                }:
                exterior = np.array([
                    [[0, np.sin(x)*r, np.cos(x)*r] for x in np.linspace(0, np.pi*2, 40)]
                    for r in np.linspace(start["t3"]/2, end["t3"]/2, 2)
                ])


class _FrameSection(_Section):
    polygon: list

    def _create(self, csi, model, library, conv):

        self.polygon = []

        section = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
                           SectionName=self.name
        )

        segments = find_rows(csi.get("FRAME SECTION PROPERTIES 05 - NONPRISMATIC",[]),
                             SectionName=section["SectionName"])

        if section is None:
            print(csi["FRAME SECTION PROPERTIES 01 - GENERAL"])
            raise Exception(f"{self.name = }")

        if section["Shape"] not in {"Nonprismatic"}:
            material = find_row(csi["MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES"],
                                Material=section["Material"]
            )

            if "G12" in material:
                model.section("FrameElastic", self.index,
                              A  = section["Area"],
                              Ay = section["AS2"],
                              Az = section["AS2"],
                              Iz = section["I33"],
                              Iy = section["I22"],
                              J  = section["TorsConst"],
                              E  = material["E1"],
                              G  = material["G12"]
                )
                self.integration.append(self.index)


        elif section["Shape"] == "Nonprismatic" and \
             len(segments) != 1: #section["NPSecType"] == "Advanced":

            # TODO: Currently just treating advanced as normal prismatic section

            assert all(segment["StartSect"] == segment["EndSect"] for segment in segments)

            if segments[0]["StartSect"] not in library: # not conv.identify("AnalSect", "section", segments[0]["StartSect"]): #
                tag = self.index # conv.define("AnalSect", "section", segments[0]["StartSect"])
                library[segments[0]["StartSect"]] = \
                        _FrameSection(segments[0]["StartSect"], csi, tag, model, library, conv)

            self.integration.append(conv.identify("AnalSect", "section", segments[0]["StartSect"]))


        # 
        elif section["Shape"] == "Nonprismatic" and \
             len(segments) == 1: #section["NPSecType"] == "Default":

            segments = find_rows(csi["FRAME SECTION PROPERTIES 05 - NONPRISMATIC"],
                                 SectionName=section["SectionName"])

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
                    start= end = material[prop]
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
            off = 1
            for x,wi in zip(*leggauss(nip)):
                xi = (1+x)/2

                model.section("FrameElastic", self.index+off, #conv.define("AnalSect", "section", None), #
                              A  = interpolate(xi, "Area"),
                              Ay = interpolate(xi, "AS2"),
                              Az = interpolate(xi, "AS2"),
                              Iz = interpolate(xi, "I33"),
                              Iy = interpolate(xi, "I22"),
                              J  = interpolate(xi, "TorsConst"),
                              E  = interpolate(xi, "E1"),
                              G  = interpolate(xi, "G12")
                )


                self.integration.append((self.index+off, xi, wi/2))

                off += 1


        else:
            # TODO: truss section?
            warnings.warn(f"Unknown shape {section['Shape']}")
            pass
