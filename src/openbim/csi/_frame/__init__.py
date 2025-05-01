#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
import sys
import numpy as np
import warnings

from ..utility import UnimplementedInstance, find_row


def _is_truss(frame, csi):
    if "FRAME RELEASE ASSIGNMENTS 1 - GENERAL" in csi:
        release = find_row(csi["FRAME RELEASE ASSIGNMENTS 1 - GENERAL"],
                        Frame=frame["Frame"])
    else:
        return False

    return release and all(release[i] for i in ("TI", "M2I", "M3I", "M2J", "M3J"))


def _orient(xi, xj, angle):
    """
    Calculate the coordinate transformation vector.
    xi is the location of node I, xj node J,
    and `angle` is the rotation about the local axis in degrees

    By default local axis 2 is always in the 1-Z plane, except if the object
    is vertical and then it is parallel to the global X axis.
    The definition of the local axes follows the right-hand rule.
    """

    # The local 1 axis points from node I to node J
    dx, dy, dz = e1 = xj - xi
    # Global z
    E3 = np.array([0, 0, 1])

    # In Sap2000, if the element is vertical, the local y-axis is the same as the
    # global x-axis, and the local z-axis can be obtained by cross-multiplying
    # the local x-axis with the local y-axis.
    if dx == 0 and dy == 0:
        e2 = np.array([1, 0, 0])

    # Otherwise, the plane composed of the local x-axis and the local
    # y-axis is a vertical plane. In this
    # case, the local z-axis can be obtained by the cross product of the local
    # x-axis and the global z-axis.
    else:
        e2 = np.cross(E3, e1)

    e3 = np.cross(e1, e2)

    # Rotate the local axis using the Rodrigue rotation formula
    # convert from degrees to radians
    angle = angle / 180 * np.pi
    e3r = e3 * np.cos(angle) + np.cross(e1, e3) * np.sin(angle)
    # Finally, the normalized local z-axis is returned
    return e3r / np.linalg.norm(e3r)


def add_frames(csi, model, library, config, conv):
    ndm = config.get("ndm", 3)
    log = []

    # itag = 1
    transform = 1

    tags = {}

    for frame in csi.get("CONNECTIVITY - FRAME",[]):
        if _is_truss(frame, csi):
            conv.log(UnimplementedInstance("Truss", frame))
            continue

        if "IsCurved" in frame and frame["IsCurved"]:
            conv.log(UnimplementedInstance("Frame.Curve", frame))

        nodes = (
            conv.identify("Joint", "node", frame["JointI"]),
            conv.identify("Joint", "node", frame["JointJ"])
        )

        if "FRAME ADDED MASS ASSIGNMENTS" in csi:
            row = find_row(csi["FRAME ADDED MASS ASSIGNMENTS"],
                            Frame=frame["Frame"])
            additional_mass = row["MassPerLen"] if row else 0.0
        else:
            additional_mass = 0.0

        #
        # Geometric transformation
        #
        if "FRAME LOCAL AXES ASSIGNMENTS 1 - TYPICAL" in csi:
            row = find_row(csi["FRAME LOCAL AXES ASSIGNMENTS 1 - TYPICAL"],
                            Frame=frame["Frame"])
            angle = row["Angle"] if row else 0.0
        else:
            angle = 0

        xi = np.array(model.nodeCoord(nodes[0]))
        xj = np.array(model.nodeCoord(nodes[1]))
        if np.linalg.norm(xj - xi) < 1e-10:
            log.append(UnimplementedInstance("Frame.ZeroLength", frame))
            print(f"ZERO LENGTH FRAME: {frame['Frame']}", file=sys.stderr)
            continue

        if ndm == 3:
            vecxz = _orient(xi, xj, angle)
            model.geomTransf("Linear", transform, *vecxz)
        else:
            model.geomTransf("Linear", transform)

        transform += 1

        #
        # Section
        #
        assign = find_row(csi["FRAME SECTION ASSIGNMENTS"], Frame=frame["Frame"])
        
        sect_info = find_row(csi["FRAME SECTION PROPERTIES 01 - GENERAL"],
                             SectionName=assign["AnalSect"])
        if not sect_info:
            conv.log(UnimplementedInstance("FrameSection.Unknown", assign))
            continue

        # ---------------------------------------------------------------------
        # Handle prismatic vs nonprismatic to get mass/length
        # ---------------------------------------------------------------------
        is_nonprismatic = (sect_info["Shape"] == "Nonprismatic")

        if not is_nonprismatic:
            #   
            # Prismatic section
            #
            A = sect_info["Area"]  # cross‐sectional area
            mat_info = find_row(csi["MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES"],
                                Material=sect_info["Material"])
            rho = mat_info["UnitMass"]  # material density (mass / volume)

            # self‐weight mass per length = area × density
            self_weight_mpl = A * rho

        else:
            #
            # Nonprismatic section
            #
            # look for total mass and total length in
            # 

            # total mass for the entire nonprismatic section
            total_mass = sect_info.get("TotalMass", 0)
            if "NPSectLen" in assign:
                total_length = assign["NPSectLen"]
            else:
                # handle the case where NPSectLen doesn't exist
                total_length = np.linalg.norm(xj - xi)

            if total_length < 1e-10:
                conv.log(UnimplementedInstance("FrameSection.NonprismaticZeroLength", assign))
                continue

            # self‐weight mass per length = total mass / total length
            self_weight_mpl = total_mass / total_length

        # Combine any assigned mass per length with self‐weight mass per length
        total_mass = self_weight_mpl + additional_mass

        # section = library["frame_sections"][assign["AnalSect"]] # conv.identify("AnalSect", "section", assign["AnalSect"]) #

        # TODO: probably dont use "SectionType" from this table; it looks like its 
        # superceded by "Shape" in Section Properties 01 table.
        if ("SectionType" not in assign) or (assign["SectionType"] != "Nonprismatic") or \
           assign["NPSectType"] == "Advanced":

            section = conv.identify("AnalSect", "section", assign["AnalSect"])
            if section is None:
                warnings.warn(f"No section found for {assign['AnalSect']}")
                continue

            e = model.element("PrismFrame",
                          conv.define("Frame", "element", frame["Frame"]),
                          nodes,
                          section=section,
                          transform=transform-1,
                          mass=total_mass
            )
            tags[frame["Frame"]] = e


        elif assign["NPSectType"] == "Default" and (integr := conv.identify("AnalSect", "integration", assign["AnalSect"])):
            # Non-prismatic sections
            e = model.element("ForceFrame",
                              conv.define("Frame", "element", frame["Frame"]),
                              nodes,
                              transform-1,
                              integr,
                              mass=total_mass
            )
            tags[frame["Frame"]] = e

        else:
            conv.log(UnimplementedInstance("FrameSection.NPSectType", assign["NPSectType"]))
            continue

    library["frame_tags"] = tags


    return log


