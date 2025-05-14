#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
import re
import sys
import json
import warnings
import numpy as np
from openbim.convert import Converter
import opensees.openseespy as ops

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
    "C3D6":  "wedge",
    "C3D15": "wedge15",
    #
    # 4-node bilinear displacement and pore pressure
    "CAX4P": "quad",
    # 6-node quadratic
    "CPE6": "triangle6",
}

meshio_to_abaqus_type = {v: k for k, v in abaqus_to_meshio_type.items()}

def _iter_nodes(block, n=None):
    if n is None:
        for line in block.data:
            if line[-1] == ",":
                line = line[:-1]
            yield map(int, re.split(",\\W*", line.strip())) # line.split(","))
    else:
        lines = iter(block.data)
        nodes = []
        while True:
            if len(nodes) == n+1:
                yield tuple(nodes)
                nodes = []
            try:
                line = next(lines)
            except StopIteration:
                if len(nodes) > 0:
                    yield tuple(nodes)
                break
            if line[-1] == ",":
                line = line[:-1]
            nodes.extend(map(int, re.split(",\\W*", line.strip())))

def _create_sections(ast, model, conv):

    for node in ast.find_all("Shell Section"):
        els = node.attributes.get("elset", None)
        tag = node.attributes.get("name")
        mat = node.attributes.get("material")
        thickness = node.attributes.get("thickness", None)

    for node in ast.find_all("Beam Section"):
        elset = node.attributes.get("elset")
        mat = node.attributes.get("material")
        shape = node.attributes.get("section", None)
        assert shape in {"RECT", "L", "T", "C", "I"}
        thickness = node.attributes.get("thickness", None)
        tag = conv.define("FrameSection", "FrameSection", elset)


def create_model(ast, verbose=False, mode=None):
    if mode is None:
        mode = "simulate"

    part = ast
    for i,node in enumerate(ast.find_all("Part")):
        if len(node.children) == 0 or i != 2:
            continue
        part = node
        break

    model, conv = _create_part(part, verbose=verbose, mode=mode)


    # Boundaries
#   _create_boundaries(ast, model, conv)

    return model


def _create_boundaries(ast, model, conv):

    for block in ast.find_all("Boundary"):
        for line in block.data:
            print(line)
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


def _create_part(ast, verbose=False, mode=None):
    if mode is None:
        mode = "simulate"

    # Create a new model
    model = ops.Model(ndm=3, ndf=6)
    conv = Converter()

    # Create Materials


    if mode == "visualize":
        # Create dummy materials to make the model work
        E = 29e3
        nu = 0.27
        density = 1.27

        mat = 1
        fsec = 1
        shell_section = 1
        model.material("ElasticIsotropic", mat, E, nu)
        model.section("ElasticShell",        1, E, nu, 1.175, density)
        model.section("FrameElastic", fsec, E=E, G=E*0.6, A=1, Iy=1, Iz=1, J=1)

        model.geomTransf("Linear", 1, (0.0, 1.0, 0))

    #
    # Parse materials
    #
    for node in ast.find_all("Material"):
        density = 0.0
        E = None
        nu = None
        Fy = None
        Ep = None
        hardening = None
        for child in node.children:
            if child.keyword == "Density":
                density = float(child.data[0].split(",")[0])
                continue

            if child.keyword == "Elastic":
                properties = child.data[0].split(",")
                E = float(properties[0])
                try:
                    nu = float(properties[1])
                except:
                    pass
                #                   model.uniaxialMaterial('Elastic', material_name, E)

            elif child.keyword == "Plastic":
                hardening = child.attributes.get("hardening", "isotropic")
                properties = child.data[0].split(",")
                if hardening == "isotropic":
                    Fy = float(properties[0])
                    Ep = float(properties[1])
                elif hardening == "JOHNSON COOK":
                    pass

            elif child.keyword == "Concrete":
                continue
                tag = conv.define("Material", "uniaxial", node.attributes.get("name"))
                properties = child.children[0].attributes.get("data").split(",")
                f_c = float(properties[0])  # Compressive strength
                f_t = float(properties[1])  # Tensile strength
                model.uniaxialMaterial("Concrete", tag, f_c, f_t)


        # Create the material
        tag = conv.define("Material", "material", node.attributes.get("name"))
        if Fy is not None:
            assert E is not None
            assert nu is not None

            if hardening == "isotropic":
                model.nDMaterial("J2Simplified", tag, 
                                 E=E, nu=nu, 
                                 Fy=Fy, 
                                 Hiso=E,
                                 Hkin=0,
                                 density=density
                )

            elif hardening == "JOHNSON COOK":
                warnings.warn("JOHNSON COOK hardening not implemented")
                model.material("ElasticIsotropic", tag, E, nu)

        elif E is not None:
            if nu is not None:
                model.material("ElasticIsotropic", tag, E, nu)
            else:
                assert E is not None
                model.uniaxialMaterial("Elastic", tag, E)

    if False:
        _create_sections(ast, model, conv)


    # Create nodes
    for nodes in ast.find_all("Node"):
        for line in nodes.data:
            node_data = line.split(",")
            node_id = int(node_data[0])
            coords = tuple(map(float, node_data[1:]))
            if len(coords) == 2:
                coords = (coords[0], coords[1], 0.0)
            elif len(coords) == 1:
                coords = (coords[0], 0.0, 0.0)
            model.node(node_id, coords)

    for nodes in ast.find_all("Ngen"):
        nset = nodes.attributes.get("nset", None)

        i_tag, j_tag = map(int, nodes.data[0].split(","))

        i_x = model.nodeCoord(i_tag)
        j_x = model.nodeCoord(j_tag)
        n  = j_tag - i_tag - 1

        for node_id, coords in zip(range(i_tag+1, j_tag), np.linspace(i_x, j_x, n)):
            coords = tuple(coords.tolist())
            if len(coords) == 2:
                coords = (coords[0], coords[1], 0.0)
            elif len(coords) == 1:
                coords = (coords[0], 0.0, 0.0)
            model.node(node_id, coords)


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

        elset = block.attributes.get("elset", None)
        section = None
        if elset is not None:
            section = None #conv.identify("FrameSection", "FrameSection", elset)

        #
        # Create the element
        #
        if element_type == "hexahedron":
            for tag, *nodes in _iter_nodes(block):
                assert len(nodes) == 8
                model.element("stdBrick", i, tuple(nodes), mat)
                i += 1

        elif element_type == "hexahedron20":
            for tag, *nodes in _iter_nodes(block, 20):
                model.element("stdBrick", i, tuple(nodes[:8]), mat)
                i += 1
                print("WARNING Brick with ", len(nodes), "nodes", file=sys.stderr)

        elif element_type == "quad":
            for tag, *nodes in _iter_nodes(block):
                if len(nodes) == 4:
                    model.element("ShellMITC4", i, list(nodes), section=shell_section)
                    i += 1
                else:
                    print("WARNING Quad with ", len(nodes), "nodes", file=sys.stderr)

        elif element_type == "quad8":
            for tag, *nodes in _iter_nodes(block):
                assert len(nodes) == 8
                model.element("Quad", i, list(nodes), section=shell_section)
                i += 1

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

        elif element_type == "line3":
            fsec = 1
            model.section("FrameElastic", fsec, E=E, G=E*0.6, A=1, Iy=1, Iz=1, J=1)
            model.geomTransf("Linear", 1, (0.0, 1.0, 0))
            for tag, *nodes in _iter_nodes(block):
                model.element("ExactFrame", i, list(nodes), section=fsec, transform=1)
                i += 1

        elif element_type == "triangle":
            for tag, *nodes in _iter_nodes(block):
                if len(nodes) == 3:
                    model.element("ShellDKGT", i, tuple(nodes), shell_section)
                    i += 1
                else:
                    print("Shell with ", len(nodes), "nodes")

        elif element_type in {"tetra", "tetra4"}:
            for tag, *nodes in _iter_nodes(block):
                assert len(nodes) == 4
                model.element("FourNodeTetrahedron", i, tuple(nodes), mat)
                i += 1

        elif element_type == "tetra10":
            for tag, *nodes in _iter_nodes(block):
                assert len(nodes) == 10
                model.element("TenNodeTetrahedron", i, tuple(nodes), mat)
                i += 1
        else:
            print("WARNING unsupported element \"", element_type, "\"", file=sys.stderr)

    return model, conv
