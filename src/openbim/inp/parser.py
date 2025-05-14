#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
# cmp
#
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
    "Part": [
        "Node", "Element", "Nset", "Elset", 
        # Sections
        "Shell Section", "Shell General Section", 
        "Beam Section", "Beam General Section",
        "Solid Section",
        "Orientation", "Surface"
    ],
    "Node": [],  # Nodes typically do not have children
    "Element": [],
    "Nset": [],
    "Elset": [],
    "Shell Section": [],
    "Beam Section": [],
    "Beam General Section": [],
    "Orientation": [],
    "Surface": [],

    "Assembly": ["Instance", "Nset", "Elset", "Surface", "Contact"],
    "Instance": [],
    "Surface": [],

    "Amplitude": [], # Amplitude seems to be a root node

    # Material properties
    "Material": [
        "Density", "Elastic", "Plastic", "Permeability", 
        "Viscoelastic", 
        "Conductivity", "Expansion",
        "Inelastic Heat Fraction",
        "Rate Dependent",
        "Specific Heat",
        "Mohr Coulomb", "Mohr Coulomb Hardening",
        "Thermal Conductivity"
    ],
    "Density": [],  # Density has no children
    "Elastic": [],  # Elastic properties have no children
    "Plastic": [],  # Plastic properties have no children
    "Conductivity": [],  # Conductivity has no children
    "Expansion": [],
    "Viscoelastic": [],  # Additional material property

    "Initial Conditions": [],

    # STEP
    "Step": [
        "Static", "Dynamic", "Heat Transfer",
        "Boundary",
        "Geostatic",
        "Field",
        "Soils",
        "Buckle",
        "Dload", "Dsload", 
        "Restart", 
        "Output", "Node Output", "Element Output", "Contact Output"
    ],

    "Static": [], 
    "Dynamic": [],  # Dynamic analysis step, no children
    "Heat Transfer": [],  # Heat transfer analysis, no children
    "Boundary": [],  # Boundary conditions, no children
    "Dload": [],  # Distributed loads, no children
    "Cload": [],  # Distributed loads, no children
    "Dsload": [],
    "Buckle": [],
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

    def add_child(self, child_node):
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

    if "GENERATE" in params_map:
        if len(set_ids) != 3:
            raise Exception(set_ids)
        set_ids = range(set_ids[0], set_ids[1] + 1, set_ids[2])
    return set_ids, set_names, line


def load(filename, verbose=False):

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
                # print(line, stack)
                n = stack.pop()
                if len(stack) == 0:
                    # print(n)
                    break
                continue

            if line.startswith("*"):  # Identify keywords
                # Split keyword and attributes
                parts = line[1:].split(",", 1)
                keyword = line.partition(",")[0].strip().replace("*", "").title() # parts[0].strip()
                attributes = {}
                if len(parts) > 1:
                    # Process attributes
                    for attr in parts[1].split(","):
                        key_value = attr.split("=")
                        if len(key_value) == 2:
                            attributes[key_value[0].strip().lower()] = key_value[1].strip()

                # Create a new node
                current_node = AbaqusTable(keyword,
                                           attributes,
                                           child_keys=hierarchy.get(keyword,[])
                )


                while len(stack) > 1:
                    if keyword in stack[-1].child_keys:
                        break
                    popped = stack.pop()
                    if verbose:
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
        if verbose:
            print(root)
        return root
    

if __name__ == "__main__":
    import sys

    ast = load(sys.argv[1])
    print(ast)
