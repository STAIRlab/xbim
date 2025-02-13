#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
import numpy as np
import warnings
from .utility import UnimplementedInstance, find_row, find_rows
from .handler import Handler


class LinkHandler(Handler):
    pass

def _orient(xi, xj, deg):
    # Calculate the direction vector of the link element
    # Where a is the node number of node i, b is the node number of node j, and degree is the user-specified local axis
    # ------------------------------------------------------------------------------
    d_x, d_y, d_z = e_x = xj - xi

    # Local 1-axis points from node I to node J
    l_x = np.array([d_x, d_y, d_z])
    # Global z-axis
    g_z = np.array([0, 0, 1])

    # In SAP2000, if the link is vertical, the local y-axis is the same as the
    # global x-axis, and the local z-axis can be obtained by crossing the local
    # x-axis with the local y-axis
    if d_x == 0 and d_y == 0:
        l_y = np.array([1, 0, 0])
        l_z = np.cross(l_x, l_y)

    # In other cases, the plane formed by the local x-axis and the local y-axis
    # is a vertical plane (i.e., the normal vector is horizontal), and the
    # local z-axis can be obtained by crossing the local x-axis with the global
    # z-axis
    else:
        l_z = np.cross(l_x, g_z)

    # The local axis may also be rotated using the Rodrigues' rotation formula
    angle = deg / 180 * np.pi
    l_z_rot = l_z * np.cos(angle) + np.cross(l_x, l_z) * np.sin(angle)
    # The rotated local y-axis can be obtained by crossing the rotated local z-axis with the local x-axis
    l_y_rot = np.cross(l_z_rot, l_x)
    # Finally, return the normalized local y-axis
    return l_y_rot / np.linalg.norm(l_y_rot)


_link_tables = {
    "Linear              " : "LINK PROPERTY DEFINITIONS 02 - LINEAR",
    "??                  " : "LINK PROPERTY DEFINITIONS 03 - MULTILINEAR",
    "Damper - Exponential" : "LINK PROPERTY DEFINITIONS 04 - DAMPER",
    "???                 " : "LINK PROPERTY DEFINITIONS 05 - GAP",
    "????                " : "LINK PROPERTY DEFINITIONS 06 - HOOK",
    "?????               " : "LINK PROPERTY DEFINITIONS 07 - RUBBER ISOLATOR",
    "??????              " : "LINK PROPERTY DEFINITIONS 08 - SLIDING ISOLATOR",
    "Plastic (Wen)"        : "LINK PROPERTY DEFINITIONS 10 - PLASTIC (WEN)",
    "???????             " : "LINK PROPERTY DEFINITIONS 11 - MULTILINEAR PLASTIC",
}

def create_links(csi, model, library, config, conv):

    for link in csi.get("CONNECTIVITY - LINK",[]):

        nodes = (
            conv.identify("Joint", "node", link["JointI"]),
            conv.identify("Joint", "node", link["JointJ"])
        )

        # if any(not (isinstance(node, int) or node.isdigit()) for node in nodes):
        #     conv.log(UnimplementedInstance(f"Joint: non-integer nodes", assign))
        #     continue

        assign = find_row(csi["LINK PROPERTY ASSIGNMENTS"],
                          Link=link["Link"])

        if assign["LinkJoints"] == "SingleJoint":

            props = find_row(csi["LINK PROPERTY DEFINITIONS 01 - GENERAL"],
                             Link=assign["LinkProp"])

            if props["LinkType"] != "Linear":
                conv.log(UnimplementedInstance(f"Joint.SingleJoint.LinkType={props['LinkType']}", assign))

            # TODO: Implement soil springs
            props = find_rows(csi["LINK PROPERTY DEFINITIONS 02 - LINEAR"],
                             Link=assign["LinkProp"])

            flags = tuple(
                1 if find_row(props, DOF=f"{dof[0]}{'XYZ'.find(dof[1])+1}") and config["dofs"][dof]
                  else 0
                  for dof in config["dofs"]
            )

            model.fix(nodes[0], flags)

            continue

        elif assign["LinkJoints"] != "TwoJoint":
            conv.log(UnimplementedInstance(f"Joint.{assign['LinkJoints']}", assign))
            continue

        #
        # Get mats and dofs
        #
        mats = tuple(library["link_materials"][assign["LinkProp"]].values())
        dofs = tuple(library["link_materials"][assign["LinkProp"]].keys())
        dofs = tuple(["U1", "U2", "U3", "R1", "R2", "R3"].index(i)+1 for i in dofs)

        if len(dofs) == 0:
            conv.log(UnimplementedInstance(f"Joint.DOFS", assign))
            continue

        # Check whether there are any zero-length link elements
        xi = np.array(model.nodeCoord(nodes[0]))
        xj = np.array(model.nodeCoord(nodes[1]))

        distance = np.linalg.norm(xj - xi)
        zero_length_threshold = 1e-6

        #
        # Get axes and orientation
        #
        axes = find_row(csi.get("LINK LOCAL AXES ASSIGNMENTS 1 - TYPICAL", []),
                        Link=link["Link"])

        orient_vector = None  # Default value

        if axes:
            if axes["AdvanceAxes"]:
                # Handle advanced axes
                axes_advance = find_row(csi.get("LINK LOCAL AXES ASSIGNMENTS 2 - ADVANCED", []),
                                        Link=link["Link"])

                # Common advanced axes setup
                orient_vector = (
                    axes_advance["AxVecX"], axes_advance["AxVecY"], axes_advance["AxVecZ"],
                    axes_advance["PlVecX"], axes_advance["PlVecY"], axes_advance["PlVecZ"],
                )

                # Additional validation for non-zero-length links
                if distance > zero_length_threshold:
                    # Calculate node-based orientation
                    orient_vector_from_nodes = _orient(xi, xj, axes["Angle"])

                    # Calculate y-axis from advanced axes
                    ax_vec = np.array([axes_advance["AxVecX"], axes_advance["AxVecY"], axes_advance["AxVecZ"]])
                    pl_vec = np.array([axes_advance["PlVecX"], axes_advance["PlVecY"], axes_advance["PlVecZ"]])

                    x_axis = ax_vec / np.linalg.norm(ax_vec)
                    pl_projection = pl_vec - np.dot(pl_vec, x_axis) * x_axis
                    y_axis_advanced = pl_projection / np.linalg.norm(pl_projection)

                    # Validate and update orientation
                    if np.allclose(y_axis_advanced, orient_vector_from_nodes, atol=1e-6):
                        orient_vector = tuple(_orient(xi, xj, axes["Angle"]))
                    else:
                        orient_vector = tuple(_orient(xi, xj, axes["Angle"]))
                        warnings.warn(f"Orientation mismatch in link {link['Link']}")
#                       raise ValueError(f"Orientation mismatch in link {link['Link']}")

            else:
                # Handle typical axes (non-advanced)
                orient_vector = tuple(_orient(xi, xj, axes["Angle"]))


        #
        # Create the link
        #
        tag = conv.define("Link", "element", link["Link"])
        if orient_vector is not None:
            model.element("TwoNodeLink",
                      tag,
                      nodes,
                      mat=mats,
                      dir=dofs,
                      orient=orient_vector
                      )
        else:
            model.element("TwoNodeLink",
                      tag,
                      nodes,
                      mat=mats,
                      dir=dofs
                      )

    return
