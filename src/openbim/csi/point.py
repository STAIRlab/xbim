#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#

from .utility import UnimplementedInstance, find_row

def create_points(csi, model, library, config, conv):
    log = []
    ndm = config["ndm"]
    ndf = config["ndf"]
#   dofs = config["dofs"]
    dofs = csi["ACTIVE DEGREES OF FREEDOM"][0]

    used = set()

    for node in csi["JOINT COORDINATES"]:

        if node["CoordSys"] != "GLOBAL":
            coordinates = tuple(node[i] if i in node else 0.0 for i in ("GlobalX", "GlobalY", "GlobalZ"))
            # log.append(UnimplementedInstance(f"Joint.Mass.CoordSys={node['CoordSys']}", node))
        else:
            coordinates = tuple(node[i] if i in node else 0.0 for i in ("XorR", "Y", "Z"))

        node_tag = conv.define("Joint", "node", node["Joint"])


        model.node(node_tag, coordinates)

        for i,v in enumerate(dofs.values()):
            if not v:
                model.fix(node_tag, dof=i+1)

    for node in csi.get("JOINT RESTRAINT ASSIGNMENTS", []):
        # # TODO!!
        # try:
        #     int(node["Joint"])
        # except:
        #     continue

        node_tag = conv.identify("Joint", "node", node["Joint"])
        # Only fix dofs that werent aready fixed globally
        # Note that dof keys look like UX, RY, etc, but in the restraint
        # table they look like U1, R2, etc
        model.fix(node_tag, tuple(int(node[f"{key[0]}{'XYZ'.find(key[1])+1}"]) if dofs[key] else 0 for key in dofs))


    for node in csi.get("JOINT ADDED MASS ASSIGNMENTS", []):
        mass = [node[f"Mass{i+1}"] for i in range(ndm)]
        mass = mass + [0.0]*(ndf-len(mass))
        node_tag = conv.identify("Joint", "node", node["Joint"])
        if node["CoordSys"] != "GLOBAL":
            conv.log(UnimplementedInstance(f"Joint.Mass.CoordSys={node['CoordSys']}", node))
        model.mass(node_tag, tuple(mass))

    
    for node in csi.get("JOINT ADDED MASS BY VOLUME ASSIGNMENTS", []):
        dens = find_row(csi.get("MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES",[]),
                        Material=node["Material"])["UnitMass"]
        vols = [node[f"Vol{i}"] for i in range(1,ndm+1) if f"Vol{i}" in node]
        vols = vols + [0.0]*(ndf-len(vols))
        mass = tuple(vol*dens for vol in vols)
        model.mass(node["Joint"], mass)

    used.add("JOINT COORDINATES")
    used.add("JOINT RESTRAINT ASSIGNMENTS")
    used.add("JOINT ADDED MASS ASSIGNMENTS")
    used.add("JOINT ADDED MASS BY VOLUME ASSIGNMENTS")


    log.extend( _apply_constraints(csi, model, library, config) )

    used.add("JOINT CONSTRAINT ASSIGNMENTS")

    return log


def _apply_constraints(sap, model, library, config):
    log = []

    # The format of body dictionary is {'node number':'constraint name'}
    constraints = {}

    for constraint in  sap.get("JOINT CONSTRAINT ASSIGNMENTS", []):
        if "Type" in constraint and constraint["Type"] == "Body":
            # map node number to constraint
            constraints[constraint["Joint"]] = constraint["Constraint"]
        else:
            log.append(UnimplementedInstance("Joint.Constraint", constraint))

    # Sort the dictionary by body name and return a list [(node, body name)]
    constraints = list(sorted(constraints.items(), key=lambda x: x[1]))


    if len(constraints) > 0:
        nodes = []
        # Assign the first body name to the pointer
        pointer = constraints[0][1]

        # Traverse the tuple. If the second element in the tuple, the body
        # name, is the same as the pointer, then store the node number, 
        # into nodes.
        for node, constraint in constraints:
            if constraint == pointer:
                nodes.append(node)
            else:
                # First write the nodes in nodes to the body file
                for le in range(len(nodes)-1):
                    model.eval(f"rigidLink beam {nodes[0]} {nodes[le + 1]}\n")
                # Restore nodes and save the node that returns False.
                nodes = []
                nodes.append(node)
                # The pointer is changed to the new body name
                pointer = constraint

        # After the for loop ends, write the nodes in the nodes of the last loop to the body file.
        for le in range(len(nodes)-1):
            model.eval(f"rigidLink beam {nodes[0]} {nodes[le + 1]}\n")

    return log
