#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#

import json
import shlex
import opensees.openseespy as ops

"""
Assembly: Contains instances of parts but does not directly contain nodes or elements.
Part: Contains nodes, elements, materials, and sections.
Material: Contains specific material definitions like elastic and plastic properties.
Section: Defines different material behaviors (e.g., elastic, plastic).
Element: Can reference sets but does not have child elements in the context of the input structure.
Step: Defines analysis steps and may include specific analysis types.
Boundary/Load: These conditions are applied to nodes or elements and don't have nested definitions.
"""


hierarchy = {
    "root": [
        "Heading", "Preprint",
        "Part", "Assembly", 
        "Material",
        "Amplitude",
        "Boundary", 
        "Initial Conditions"
    ],


    "Heading": [],  # General information, no children
    "Preprint": [],  # Print control, no children
    "Part": ["Node", "Element", "Nset", "Elset", 
             "Shell Section", "Beam Section", "Solid Section"
    ],
    "Node": [],  # Nodes typically do not have children
    "Element": [],
    "Nset": [],
    "Elset": [],
    "Shell Section": [],
    "Beam Section": [],
    "Beam General Section": [],

    "Assembly": ["Instance", "Nset", "Elset", "Surface", "Contact"],
    "Instance": [],
    "Surface": [],

    "Amplitude": [], # Amplitude seems to be a root node

    # Material properties
    "Material": ["Density", "Elastic", "Plastic", "Permeability", 
                 "Viscoelastic", 
                 "Mohr Coulomb", "Mohr Coulomb Hardening",
                 "Thermal Conductivity"],
    "Density": [],  # Density has no children
    "Elastic": [],  # Elastic properties have no children
    "Plastic": [],  # Plastic properties have no children
    "Viscoelastic": [],  # Additional material property

    "Initial Conditions": [],

    # STEP
    "Step": ["Static", "Dynamic", "Heat Transfer",
             "Boundary",
             "Geostatic",
             "Field",
             "Soils",
             "Dload", "Dsload", 
             "Restart", "Output"],

    "Static": [],  # Static analysis step, no children
    "Dynamic": [],  # Dynamic analysis step, no children
    "Heat Transfer": [],  # Heat transfer analysis, no children
    "Boundary": [],  # Boundary conditions, no children
    "Dload": [],  # Distributed loads, no children
    "Cload": [],  # Distributed loads, no children
    "Dsload": [],
    "Restart": [],  # Restart options, no children

    "Output": [],  # Output requests
#   "Field Output": [],  # Field output requests, no children
#   "History Output": [],
    "Load": [],  # General load definitions, no children

    "Contact": ["Interaction"],  # Contact definitions
    "Interaction": [],  # Interaction definitions, no children
}
hierarchy.update({"root": list(hierarchy.keys())})
hierarchy["Instance"] = hierarchy["Part"]


class AbaqusTable:
    def __init__(self, keyword: str, attributes: dict = None, child_keys=None):
        self.keyword = keyword
        self.attributes = attributes if attributes else {}
        self.children = []
        self.data = []
        self.child_keys = child_keys

    def add_child(self, child_node: "Node"):
        self.children.append(child_node)

    def find_attr(self, keyword, **attrs):
        for node in self.find_all(keyword):
            for attr in attrs:
                if attr not in node.attributes or (
                        node.attributes[attr] != attrs[attr]):
                    break

            else:
                return node


    def find_all(self, keyword):
        for child in self.children:
            if child.keyword == keyword:
                yield child
            else:
                yield from child.find_all(keyword)

    def _open_tag(self):
        tag = f"<{self.keyword}"

        if "name" in self.attributes:
            tag += f" name={self.attributes['name']}"

        if len(self.children):
            tag += ">\n"

        else:
            tag += " />\n"

        return tag

    def __repr__(self, level=0):
        # return f"<{self.keyword}>"
        ret = "  " * level + self._open_tag()
        for child in self.children:
            ret += child.__repr__(level + 1)

        if len(self.children):
            ret += "  " * level + f"</{self.keyword}>\n"

        return ret

def _read_set(f, params_map):
    """
    From meshio
    """
    set_ids = []
    set_names = []
    while True:
        line = f.readline()
        if not line or line.startswith("*"):
            break
        if line.strip() == "":
            continue

        line = line.strip().strip(",").split(",")
        if line[0].isnumeric():
            set_ids += [int(k) for k in line]
        else:
            set_names.append(line[0])

    set_ids = np.array(set_ids, dtype="int32")
    if "GENERATE" in params_map:
        if len(set_ids) != 3:
            raise ReadError(set_ids)
        set_ids = np.arange(set_ids[0], set_ids[1] + 1, set_ids[2], dtype="int32")
    return set_ids, set_names, line


def load(filename):

    with open(filename, "r") as file:
        root = current_node = AbaqusTable("root", child_keys=hierarchy["root"])
        stack = [current_node]

        for line in file:
            line = line.strip()
            if not line:
                # Skip empty lines
                continue

            if line.startswith("#") or line.startswith("**"):
                # Skip comments
                continue

            if line.startswith("*End"):
                stack.pop()
                continue

            if line.startswith("*"):  # Identify keywords
                # Split keyword and attributes
                parts = line[1:].split(",", 1)
                keyword = line.partition(",")[0].strip().replace("*", "").title() # parts[0].strip()
                print(keyword)
                attributes = {}
                if len(parts) > 1:
                    # Process attributes
                    for attr in parts[1].split(","):
                        key_value = attr.split("=")
                        if len(key_value) == 2:
                            attributes[key_value[0].strip()] = key_value[1].strip()

                # Create a new node
                current_node = AbaqusTable(keyword,
                                           attributes,
                                           child_keys=hierarchy.get(keyword,[])
                )


                while len(stack) > 1:
                    if keyword in stack[-1].child_keys:
                        break
                    popped = stack.pop()
                    if True:
                        print(f">> Popped {popped.keyword} from {keyword}; parent is {stack[-1].keyword}")
                stack[-1].add_child(current_node)

#               # Add to the parent node
#               if stack[-1].child_keys and keyword in stack[-1].child_keys:
#                   stack[-1].add_child(current_node)

#               elif len(stack) > 1:
#                   # Close the current parent
#                   popped = stack.pop()
#                   if True:
#                       print(f">> Popped {popped.keyword} from {keyword}; parent is {stack[-1].keyword}")


                # Check if this keyword has children
#               # Add to stack if has children
                if current_node.child_keys:
                    stack.append(current_node)


            elif current_node:
                current_node.data.append(line)

        return root


def create_opensees_model(ast):
    # Create a new model
    model = ops.Model(ndm=3, ndf=3)

    # Dictionary to map material names/IDs
    material_map = {}
    section_map = {}

    # Parse materials
#   for node in ast.find_all("Material"):
#       for child in node.children:
#           if child.keyword == "Elastic":
#               material_name = node.attributes.get("name")
#               properties = child.children[0].attributes.get("data").split(",")
#               E = float(properties[0])
#               nu = float(properties[1])
#               #                   model.uniaxialMaterial('Elastic', material_name, E)
#               model.material("ElasticIsotropic", material_name, E, nu)

#           elif child.keyword == "Plastic":
#               material_name = node.attributes.get("name")
#               properties = child.children[0].attributes.get("data").split(",")
#               E = float(properties[0])
#               yield_strength = float(properties[1])
#               model.uniaxialMaterial("Plastic", material_name, E, yield_strength)

#           elif child.keyword == "Concrete":
#               material_name = node.attributes.get("name")
#               properties = child.children[0].attributes.get("data").split(",")
#               f_c = float(properties[0])  # Compressive strength
#               f_t = float(properties[1])  # Tensile strength
#               model.uniaxialMaterial("Concrete", material_name, f_c, f_t)

#           material_map[node.attributes.get("name")] = material_name

#       if node.keyword == "Section":
#           section_name = node.attributes.get("name")
#           material_name = node.attributes.get("material")
#           thickness = node.attributes.get(
#               "thickness", None
#           )  # Optional, for 2D elements

#           # Store the section information
#           section_map[section_name] = {
#               "material": material_name,
#               "thickness": thickness,
#           }

    if child.keyword == "Node":
        for entry in child.children:
            node_data = entry.attributes.get("data").split(",")
            node_id = int(node_data[0])
            coords = list(map(float, node_data[1:]))
            model.node(node_id, *coords)

    # Create elements
    for node in ast.children:
        if node.keyword == "Part":
            for child in node.children:

                if child.keyword == "Element":
                    for entry in child.children:
                        element_data = entry.attributes.get("data").split(",")
                        element_id = int(element_data[0])
                        element_type = element_data[1]
                        connectivity = list(map(int, element_data[2:]))

                        # Extract material assignment
                        material_id = entry.attributes.get(
                            "material"
                        )  # Assume material ID is stored here

                        # Tetrahedral elements
                        if element_type == "C3D4":
                            model.element("tetrahedron", element_id, *connectivity)

                        # BEAMS
                        elif element_type == "B31":
                            model.element(
                                "elasticBeamColumn", element_id, *connectivity
                            )  # Linear 2-node beam

                        elif element_type == "B32":
                            model.element(
                                "elasticBeamColumn", element_id, *connectivity
                            )  # Quadratic 3-node beam

                        elif element_type == "B33":
                            model.element(
                                "elasticBeamColumn", element_id, *connectivity
                            )  # Linear 3-node beam

                        elif element_type == "B21":
                            model.element(
                                "elasticBeamColumn", element_id, *connectivity
                            )  # 2D beam

                        elif element_type == "B22":
                            model.element(
                                "elasticBeamColumn", element_id, *connectivity
                            )  # Quadratic 2-node beam

                        elif element_type == "C0D2":  # Beam example
                            model.element(
                                "elasticBeamColumn", element_id, *connectivity
                            )

                        # SOLID
                        elif element_type == "C3D8":  # Hexahedral element
                            model.element("brick", element_id, *connectivity)

                        elif element_type == "C2D4":  # 2D quadrilateral element
                            model.element("quad", element_id, *connectivity)

                        elif (
                            element_type == "C3D10"
                        ):  # Tetrahedral element with mid-side nodes
                            model.element(
                                "tetrahedron", element_id, *connectivity
                            )  # or use specific type

                        else:
                            print(
                                f"Warning: Unrecognized element type {element_type} for element ID {element_id}."
                            )


        elif node.keyword == "Load":
            for child in node.children:
                load_data = child.attributes.get("data").split(",")
                node_id = int(load_data[0])
                load_values = list(map(float, load_data[1:]))
                model.load(node_id, *load_values)


if __name__ == "__main__":
    import sys

    ast = load(sys.argv[1])
    print(ast)

#   create_opensees_model(ast)
