#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
"""
"SECTION DESIGNER PROPERTIES 01 - GENERAL"
"SECTION DESIGNER PROPERTIES 02 - REINFORCING AT SHAPE EDGES"
"SECTION DESIGNER PROPERTIES 03 - REINFORCING AT SHAPE CORNERS"
"SECTION DESIGNER PROPERTIES 04 - SHAPE I/WIDE FLANGE"
"SECTION DESIGNER PROPERTIES 05 - SHAPE CHANNEL"
"SECTION DESIGNER PROPERTIES 06 - SHAPE TEE"
"SECTION DESIGNER PROPERTIES 07 - SHAPE ANGLE"
"SECTION DESIGNER PROPERTIES 08 - SHAPE DOUBLE ANGLE"
"SECTION DESIGNER PROPERTIES 09 - SHAPE BOX/TUBE"
"SECTION DESIGNER PROPERTIES 10 - SHAPE PIPE"
"SECTION DESIGNER PROPERTIES 11 - SHAPE PLATE"
"SECTION DESIGNER PROPERTIES 12 - SHAPE SOLID RECTANGLE"
"SECTION DESIGNER PROPERTIES 13 - SHAPE SOLID CIRCLE"
"SECTION DESIGNER PROPERTIES 14 - SHAPE SOLID SEGMENT"
"SECTION DESIGNER PROPERTIES 15 - SHAPE SOLID SECTOR"
"SECTION DESIGNER PROPERTIES 16 - SHAPE POLYGON"
"SECTION DESIGNER PROPERTIES 17 - SHAPE REINFORCING SINGLE"
"SECTION DESIGNER PROPERTIES 18 - SHAPE REINFORCING LINE"
"SECTION DESIGNER PROPERTIES 19 - SHAPE REINFORCING RECTANGLE"
"SECTION DESIGNER PROPERTIES 20 - SHAPE REINFORCING CIRCLE"
"SECTION DESIGNER PROPERTIES 21 - SHAPE REFERENCE LINE"
"SECTION DESIGNER PROPERTIES 22 - SHAPE REFERENCE CIRCLE"
"SECTION DESIGNER PROPERTIES 23 - SHAPE CALTRANS SQUARE"
"SECTION DESIGNER PROPERTIES 24 - SHAPE CALTRANS CIRCLE"
"SECTION DESIGNER PROPERTIES 25 - SHAPE CALTRANS HEXAGON"
"SECTION DESIGNER PROPERTIES 26 - SHAPE CALTRANS OCTAGON"
"SECTION DESIGNER PROPERTIES 27 - CALTRANS LONGITUDINAL REBAR"
"SECTION DESIGNER PROPERTIES 28 - CALTRANS LONGITUDINAL TENDONS"
"SECTION DESIGNER PROPERTIES 29 - CALTRANS CONFINEMENT REBAR"
"SECTION DESIGNER PROPERTIES 30 - FIBER GENERAL"
"SECTION DESIGNER PROPERTIES 31 - STRESS POINT"
"SECTION DESIGNER PROPERTIES 32 - BRIDGE SECTION SHELL LAYOUT"
"SECTION DESIGNER PROPERTIES 33 - BRIDGE SECTION SOLID LAYOUT"
"SECTION DESIGNER PROPERTIES 34 - BRIDGE SECTION CUTLINE"
"SECTION DESIGNER PROPERTIES 35 - BRIDGE SECTION CENTERLINE"
"""
from ..utility import find_row, find_rows, UnimplementedInstance
import numpy as np
import warnings
from ...convert import Converter
from xbim.frame.section import ElasticConstants, PlasticConstants, FrameSection
from xbim.frame.shapes import Rectangle, WideFlange
from xsection import PolygonSection as Polygon

CIRCLE_DIVS = 40


def add_frame_sections(csi, model, conv):
    for sect in csi.get("FRAME SECTION PROPERTIES 01 - GENERAL", []):

        if not conv.identify("AnalSect", "section",     sect["SectionName"]) and \
           not conv.identify("AnalSect", "integration", sect["SectionName"]):

            if (s:= create_section(csi, sect, conv)) is not None:
                if (e := s.elastic()) is not None:
                    model.section("FrameElastic",
                                conv.define("AnalSect", "section", sect["SectionName"]), #self.index,
                                A  = e.A,
                                Ay = e.Ay,
                                Az = e.Az,
                                Iz = e.Iz,
                                Iy = e.Iy,
                                J  = e.J,
                                E  = e.E,
                                G  = e.G
                    )
                continue

            if (sections := _create_integration(csi, sect)) is not None:
                tags = []
                for s in sections:
                    tag = conv.define("AnalSect", "section", sect["SectionName"])
                    tags.append(tag)

                    if (e := s[0].elastic()) is not None:
                        model.section("FrameElastic",
                                    tag,
                                    A  = e.A,
                                    Ay = e.Ay,
                                    Az = e.Az,
                                    Iz = e.Iz,
                                    Iy = e.Iy,
                                    J  = e.J,
                                    E  = e.E,
                                    G  = e.G
                        )

                model.beamIntegration("UserDefined",
                          conv.define("AnalSect", "integration", sect["SectionName"]),
                          len(sections),
                          tuple(tags),
                          tuple(i[1] for i in sections),
                          tuple(i[2] for i in sections))
                continue

            conv.log(UnimplementedInstance(f"FrameSection.Shape={sect['Shape']}"))
            # assert False, sect

    return


def iter_sections(csi, conv=None):
    if conv is None:
        conv = Converter()

    for sect in csi.get("FRAME SECTION PROPERTIES 01 - GENERAL", []):

        if not conv.identify("AnalSect", "section",     sect["SectionName"]) and \
           not conv.identify("AnalSect", "integration", sect["SectionName"]):

            if (s:= create_section(csi, sect, conv)) is not None:
                yield s, sect["SectionName"]

            else:
                conv.log(UnimplementedInstance(f"FrameSection.Shape={sect['Shape']}"))
                continue
                assert False, sect



def create_section(csi, prop_01, conv=None,
                 elastic_only=False,
                 render_only=False) ->"_CsiSection":
    if conv is None:
        conv = Converter()

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

    #
    # 1)
    #

    segments = find_rows(csi.get("FRAME SECTION PROPERTIES 05 - NONPRISMATIC",[]),
                            SectionName=prop_01["SectionName"])


    if prop_01["Shape"] not in {"Nonprismatic"}:
        s = _CsiSection(csi, prop_01)
        return s

    # Note: if len(segments) == 1 it will be handled in _create_integration
    # 2)
    if prop_01["Shape"] == "Nonprismatic" and len(segments) != 1 and all(segment["StartSect"] == segment["EndSect"] for segment in segments): #section["NPSecType"] == "Advanced":

        # TODO: Currently just treating advanced as normal prismatic section

        if not conv.identify("AnalSect", "section", prop_01["SectionName"]) : #segments[0]["StartSect"]):


            # find properties
            p = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
                                SectionName=segments[0]["StartSect"])

            assert p is not None

            return _CsiSection(csi, p)


def _create_integration(csi, prop_01):
    # 3)
    segments = find_rows(csi["FRAME SECTION PROPERTIES 05 - NONPRISMATIC"],
                            SectionName=prop_01["SectionName"])

    if prop_01["Shape"] != "Nonprismatic" or len(segments) != 1: 
        return None


    # Interpolate linear-elastic sections properties
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
        assert material == find_row(csi["MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES"],
                            Material=sj["Material"]
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

    #
    # Define a numerical integration scheme
    #
    from numpy.polynomial.legendre import leggauss
    nip = 5
    sections = []
    for x,wi in zip(*leggauss(nip)):
        xi = (1+x)/2

        sec = _InterpolatedSection(
            ElasticConstants(
                        A  = interpolate(xi, "Area"),
                        Ay = interpolate(xi, "AS2"),
                        Az = interpolate(xi, "AS2"),
                        Iz = interpolate(xi, "I33"),
                        Iy = interpolate(xi, "I22"),
                        J  = interpolate(xi, "TorsConst"),
                        E  = interpolate(xi, "E1"),
                        G  = interpolate(xi, "G12")
            )
        )

        sections.append((sec, xi, wi/2))


    return sections



class _InterpolatedSection: #(FrameSection):
    def __init__(self, elastic):
        self._elastic = elastic
        self._plastic = None

    def elastic(self):
        return self._elastic

    def plastic(self):
        return self._plastic

    def create_fibers(self, mesh_scale=None, **kwds):
        return None

class _CsiSection: # (FrameSection):
    def __init__(self, csi, prop_01, mode=None):
        self._mode = mode

        if isinstance(prop_01, str):
            name = prop_01
            prop_01 = find_row(csi.get("FRAME SECTION PROPERTIES 01 - GENERAL",[]), SectionName=name)
            if prop_01 is None:
                raise ValueError(f"Section {name} not found.")
        else:
            name = prop_01["SectionName"]

        self._prop_01 = prop_01
        self._csi = csi
        self._offset = prop_01.get("CGOffset2", 0), prop_01.get("CGOffset3", 0)

        self._elastic = None

    # Simulation
    def create_fibers(self, mesh_scale=None, **kwds):
        return self._create_model(mesh_scale=mesh_scale).create_fibers(**kwds)

    def plastic(self):
        pass

    def elastic(self):
        prop_01 = self._prop_01
        csi = self._csi

        if self._elastic is None:
            self._elastic = _find_elastic(csi, prop_01)

        return self._elastic


    def _create_model(self, mesh_size=None, render_only=True):
        csi = self._csi
        prop_01 = self._prop_01
        name = self._prop_01["SectionName"]
        details = not render_only

        centroid = [prop_01.get("CGOffset2",0), prop_01.get("CGOffset3", 0)]

        #
        # First skip over non-basic sections if we need detailed model
        #

        shape = self._prop_01.get("Shape", "")

        if details and (shape in {"SD Section", "Bridge Section"}):
            pass

        elif details and find_row(csi.get("FRAME SECTION PROPERTIES 02 - CONCRETE COLUMN",[]), SectionName=name):
            pass


        #
        # Begin matching homogeneous standard shapes
        #
        elif shape == "Circle":
            r = prop_01["t3"]/2
            exterior = np.array([
                [np.sin(x)*r, np.cos(x)*r] for x in np.linspace(0, np.pi*2, CIRCLE_DIVS)
            ])
            return Polygon(exterior, #elastic=self.elastic(), 
                           mesh_size=r/5)

        elif shape == "I/Wide Flange":
            # Example 1-018a.s2k

            return WideFlange(
                d=prop_01["t3"],
                b=prop_01["t2"],
                tf=prop_01["tf"],
                tw=prop_01["tw"],
                centroid=centroid,
                #elastic = self.elastic()
            )
        

        elif shape == "Rectangular":
            from xbim.frame.shapes import Rectangle
            centroid = [prop_01.get("CGOffset2",0), prop_01.get("CGOffset3", 0)]
            return Rectangle(d=prop_01["t3"], 
                             b=prop_01["t2"],
                             centroid=centroid,
                             #elastic = self.elastic()
            )
        
        elif shape=="Box/Tube":
            from xbim.frame.shapes import HollowRectangle

            "SECTION DESIGNER PROPERTIES 09 - SHAPE BOX/TUBE"
            # t3=4.5   t2=4.5   tf=0.8   tw=0.8
            return HollowRectangle(
                d=prop_01["t3"],
                b=prop_01["t2"],
                tf=prop_01["tf"],
                tw=prop_01["tw"],
                centroid=centroid,
                # elastic = self.elastic()
            )
        

        elif shape == "SD Section":
            prop_sd = find_row(csi.get("SECTION DESIGNER PROPERTIES 01 - GENERAL", []), SectionName=name)

            if prop_sd["nCaltransCr"] == 1:
                circle = find_row(csi.get("SECTION DESIGNER PROPERTIES 24 - SHAPE CALTRANS CIRCLE", []), SectionName=name)
                if details:
                    return None
    
                assert circle["Height"] == circle["Width"]
                r = circle["Height"]/2
                exterior = np.array([
                    [np.sin(x)*r, np.cos(x)*r] for x in np.linspace(0, np.pi*2, CIRCLE_DIVS)
                ])
                return Polygon(exterior, 
                               mesh_size=r/5,
                            #    centroid=centroid, 
                            #    elastic=self.elastic()
                )


            elif prop_sd["nPolygon"] > 0:
                # if mesh_size is None:
                #     raise ValueError("Mesh size must be provided for polygonal sections")
                
                interior = []
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

                return Polygon(exterior, interior, 
                               mesh_size=mesh_size,
                            #    centroid=centroid, elastic=self.elastic()
                               )


        #
        # BRIDGE SECTIONS
        #
        elif shape == "Bridge Section":

            polygon_data = csi.get("FRAME SECTION PROPERTIES 06 - POLYGON DATA", [])

            cbg = find_row(csi.get("BRIDGE SECTION DEFINITIONS 02 - CONCRETE BOX GIRDER", []), Section=name)
            if cbg:
                """
                Section=BSEC1   
                ExtGirdType=AASHTO-PCI-ASBI  
                    TotalWidth=330.708661417323
                    TotalDepth=70.8661417322835
                    ExtGRadClip=22.8346456692913   
                    GirdVert=No   
                    TopSlabThk=8.85826771653543   BotSlabThk=7.8740157480315
                    ExtGirdThk=10.5905511811024   TopSlabProf=Flat   
                    LOLength=53.9370078740157   LOOuterThk=8.85826771653543   ROLength=53.9370078740157   ROOuterThk=8.85826771653543   
                    FHoriz1=53.9370078740157   FHoriz2=47.244094488189   FHoriz3=0
                    FVert1=4.92125984251969   FVert2=4.92125984251969   FVert3=0   
                    CurbLoc="Program Determined"   CurbLeft=15   CurbRight=15   CurbMedian=0   CurbMWidth=0   RefPtOffX=0   RefPtOffY=0   
                    DsgnTopDis=13.7795275590551   DsgnBotDis=7.8740157480315
                """

            elif polygon_data and  not details:
                if mesh_size is None:
                    if render_only:
                        mesh_size=10
                    else:
                        raise ValueError("Mesh size must be provided for polygonal sections")
                interior = []
                exterior_row = find_row(polygon_data, SectionName = name, Opening=False)
                exterior =  np.array([
                    [row["X"], row["Y"]]
                    for row in find_rows(polygon_data, SectionName = name) if row["Polygon"] == exterior_row["Polygon"]
                ])
                ref = (exterior_row["RefPtX"],  exterior_row["RefPtY"])

                for i in range(len(exterior)):
                    exterior[i] -= ref


                for hole in find_rows(polygon_data, SectionName = name, Opening=True):
                    interior.append(np.array([
                        [row["X"], row["Y"]]
                        for row in find_rows(polygon_data, Polygon=hole["Polygon"])
                    ]))
                    for i in range(len(interior[-1])):
                        interior[-1][i] -= ref

                return Polygon(exterior, interior, 
                               mesh_size=mesh_size,
                            #    centroid=centroid,
                            #    elastic=self.elastic()
                               )

#
#
#
def _find_elastic(csi, prop_01):
    material = find_row(csi.get("MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES", []),
                        Material=prop_01["Material"]
    )

    # assert material is not None, f"Material {prop_01['Material']} not found in MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES"
    if material is None:
        print(prop_01)

    if "G12" in material:
        return ElasticConstants(
                            A  = prop_01["Area"],
                            Ay = prop_01["AS3"],
                            Az = prop_01["AS2"], # TODO: check this!!!
                            Iy = prop_01["I22"],
                            Iz = prop_01["I33"],
                            J  = prop_01["TorsConst"],
                            E  = material["E1"],
                            G  = material["G12"]
            )
