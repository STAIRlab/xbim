#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
import opensees.openseespy as ops
import json
import sys
from openbim.convert import Converter

abaqus_to_meshio_type = {
    # trusses
    "T2D2":  "line",
    "T2D2H": "line",
    "T2D3":  "line3",
    "T2D3H": "line3",
    "T3D2":  "line",
    "T3D2H": "line",
    "T3D3":  "line3",
    "T3D3H": "line3",
    # beams
    "B21":   "line",
    "B21H":  "line",
    "B22":   "line3",
    "B22H":  "line3",
    "B31":   "line",
    "B31H":  "line",
    "B32":   "line3",
    "B32H":  "line3",
    "B33":   "line3",
    "B33H":  "line3",
    # surfaces
    "CPS4":  "quad",
    "CPS4R": "quad",
    "S4":    "quad",
    "S4R":   "quad",
    "S4RS":  "quad",
    "S4RSW": "quad",
    "S4R5":  "quad",
    "S8R":   "quad8",
    "S8R5":  "quad8",
    "S9R5":  "quad9",
    #
    "CPS3":  "triangle",
    "STRI3": "triangle",
    "S3":    "triangle",
    "S3R":   "triangle",
    "S3RS":  "triangle",
    "R3D3":  "triangle",
    #
    "STRI65": "triangle6",
    # 'TRISHELL6': 'triangle6',
    # volumes
    "C3D8":   "hexahedron",
    "C3D8H":  "hexahedron",
    "C3D8I":  "hexahedron",
    "C3D8IH": "hexahedron",
    "C3D8R":  "hexahedron",
    "C3D8RH": "hexahedron",
    # "HEX9": "hexahedron9",
    "C3D20":  "hexahedron20",
    "C3D20H": "hexahedron20",
    "C3D20R": "hexahedron20",
    "C3D20RH": "hexahedron20",
    # "HEX27": "hexahedron27",
    #
    "C3D4": "tetra",
    "C3D4H": "tetra4",
    # "TETRA8": "tetra8",
    "C3D10": "tetra10",
    "C3D10H": "tetra10",
    "C3D10I": "tetra10",
    "C3D10M": "tetra10",
    "C3D10MH": "tetra10",
    # "TETRA14": "tetra14",
    #
    # "PYRAMID": "pyramid",
    "C3D6": "wedge",
    "C3D15": "wedge15",
    #
    # 4-node bilinear displacement and pore pressure
    "CAX4P": "quad",
    # 6-node quadratic
    "CPE6": "triangle6",
}

meshio_to_abaqus_type = {v: k for k, v in abaqus_to_meshio_type.items()}

def _iter_nodes(block):
    for line in block.data:
        yield map(int, line.split(","))

def create_model(ast, verbose=False):
    # Create a new model
    model = ops.Model(ndm=3, ndf=6)
    conv = Converter()

    # Create Materials

    E = 29e3
    nu = 0.2
    mat = 1
    fsec = 1
    model.material("ElasticIsotropic", mat, E, nu)
    #                                       secTag  E     nu     h    rho
    model.section("ElasticMembranePlateSection", 1, E, 0.25, 1.175, 1.27)

    model.section("FrameElastic", fsec, E=E, G=E*0.6, A=1, Iy=1, Iz=1, J=1)

    model.geomTransf("Linear", 1, (0.0, 1.0, 0))

    #

    # Parse materials
    if True:
        # Dictionary to map material names/IDs
        section_map = {}

        for node in ast.find_all("Material"):
            for child in node.children:
                if child.keyword == "Elastic":
                    tag = conv.define("Material", "material", node.attributes.get("name"))
                    properties = child.data[0].split(",")
                    E = float(properties[0])
                    nu = float(properties[1])
                    #                   model.uniaxialMaterial('Elastic', material_name, E)
                    model.material("ElasticIsotropic", tag, E, nu)

                elif child.keyword == "Plastic":
                    continue
                    tag = conv.define("Material", "material", node.attributes.get("name"))
                    properties = child.data[0].split(",")
                    E = float(properties[0])
                    yield_strength = float(properties[1])
                    model.uniaxialMaterial("Plastic", tag, E, yield_strength)

                elif child.keyword == "Concrete":
                    continue
                    tag = conv.define("Material", "uniaxial", node.attributes.get("name"))
                    properties = child.children[0].attributes.get("data").split(",")
                    f_c = float(properties[0])  # Compressive strength
                    f_t = float(properties[1])  # Tensile strength
                    model.uniaxialMaterial("Concrete", tag, f_c, f_t)


            if node.keyword == "Section":
                tag = node.attributes.get("name")
                tag = node.attributes.get("material")
                thickness = node.attributes.get("thickness", None)

                # Store the section information
                section_map[tag] = {
                    "material": tag,
                    "thickness": thickness,
                }

    for nodes in ast.find_all("Node"):
        for line in nodes.data:
            node_data = line.split(",")
            node_id = int(node_data[0])
            coords = tuple(map(float, node_data[1:]))
            model.node(node_id, coords)


    for block in ast.find_all("Boundary"):
        for line in block.data:
            try:
                boundary_data = json.loads("["+line+"]") # line.split(",")
                dofs = tuple(map(int, boundary_data[1:]))
            except:
                print("WARNING ", line, file=sys.stderr)
                continue
            if len(dofs) > 1:
                dofs = tuple(range(dofs[0], dofs[1]+1))

            try:
                nodes = (int(boundary_data[0]), )
            except:
                # its a set name
                nodes = (
                        int(i)
                        for row in ast.find_attr("Nset", nset=boundary_data[0]).data
                        for i in row.split(",")
                )
            nodes = list(nodes)

            for node in nodes:
                for dof in dofs:
                    model.fix(node, dof=dof)

    if False:

        if node.keyword == "Load":
            for child in node.children:
                load_data = child.attributes.get("data").split(",")
                node_id = int(load_data[0])
                load_values = list(map(float, load_data[1:]))
                model.load(node_id, *load_values)

    #
    # Create elements
    #
    i = 1
    for block in ast.find_all("Element"):
        try:
            element_type = abaqus_to_meshio_type[block.attributes.get("type", block.attributes.get("TYPE", None))]
        except Exception as e:
            print("WARNING ", block.attributes, e, file=sys.stderr)
            continue

        if element_type == "hexahedron":
            for tag, *nodes in _iter_nodes(block):
                if len(nodes) == 8:
                    model.element("stdBrick", i, tuple(nodes), mat)
                    i += 1
                else:
                    print("WARNING Brick with ", len(nodes), "nodes", file=sys.stderr)

        elif element_type == "quad":
            for tag, *nodes in _iter_nodes(block):
                if len(nodes) == 4:
                    model.element("ShellMITC4", i, tuple(nodes), mat, 12.0)
                    i += 1
                else:
                    print("WARNING Quad with ", len(nodes), "nodes", file=sys.stderr)

        elif element_type == "line":
            fsec = 1
            model.section("FrameElastic", fsec, E=E, G=E*0.6, A=1, Iy=1, Iz=1, J=1)
            model.geomTransf("Linear", 1, (0.0, 1.0, 0))
            for tag, *nodes in _iter_nodes(block):
                if len(nodes) == 2:
                    model.element("PrismFrame", i, tuple(nodes), section=fsec, transform=1)
                    i += 1
                else:
                    print("Frame with ", len(nodes), "nodes", file=sys.stderr)

        elif element_type == "triangle":
            for tag, *nodes in _iter_nodes(block):
                if len(nodes) == 3:
                    model.element("ShellDKGT", i, tuple(nodes), mat, 12.0)
                    i += 1
                else:
                    print("Shell with ", len(nodes), "nodes")
        else:
            print(element_type)

        if False:
            for entry in block.data:
                element_data = entry.split(",")
                element_id   = int(element_data[0])
                nodes = tuple(map(int, element_data[1:]))

                # TODO: Extract material assignment
                mat = 1
                sec = 1
                trn = 1

                # Tetrahedral elements
                if element_type == "C3D4":
                    model.element("tetrahedron", element_id, nodes)

                # BEAMS
                elif element_type == "B31":
                    model.element("PrismFrame", element_id, nodes)  # Linear 2-node beam

                elif element_type == "B32":
                    model.element("PrismFrame", element_id, nodes)  # Quadratic 3-node beam

                elif element_type == "B33":
                    model.element("PrismFrame", element_id, nodes)  # Linear 3-node beam

                elif element_type == "B21":
                    model.element("PrismFrame", element_id, nodes)  # 2D beam

                elif element_type == "B22":
                    # Quadratic 2-node beam
                    model.element("PrismFrame", element_id, nodes, section=sec, transform=trn)


                # SOLID
                elif element_type == "C3D8":  # Hexahedral element
                    model.element("brick", element_id, nodes, mat)

                elif element_type == "C2D4":  # 2D quadrilateral element
                    model.element("quad", element_id, nodes)

                elif element_type == "C3D10":
                    # Tetrahedral element with mid-side nodes
                    model.element("tetrahedron", element_id, nodes)

                else:
                    print(
                        f"Warning: Unrecognized element type {element_type} for element ID {element_id}."
                    )
    return model
