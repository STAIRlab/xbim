#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
# Certain operations are loosley adapted from:
#    https://github.com/XunXun-Zhou/Sap2OpenSees/blob/main/STO_ver1.0.py
#
#
import numpy as np
from ..convert import Converter
from .parse import load
from .utility import UnimplementedInstance, print_log
from ._frame import create_frames
from ._shell import create_shells
from .point import create_points
from .link import create_links
from ._section import (
    create_shell_sections
)
from ._frame.section import (
    create_frame_sections, 
    section_geometry,
    collect_geometry as collect_outlines
)

CONFIG = {
    "Frame": {
        "Taper": "Subdivide", # Integrate
        "Element": "PrismFrame",
    }
}


def create_materials(csi, model, conv):
    library = conv._library

    # 1) Material

    #
    # 2) Links
    #
    mat_total = 1

    for link in csi.get("LINK PROPERTY DEFINITIONS 02 - LINEAR", []):
        if link["Fixed"]:
            conv.log(UnimplementedInstance("Link.Fixed", link))
            pass

        name = link["Link"]
        if "R" in link["DOF"]:
            stiff = link["RotKE"]
            damp  = link["RotCE"]
        else:
            stiff = link["TransKE"]
            damp  = link["TransCE"]

        # TODO: use damp
        model.eval(f"material ElasticIsotropic {mat_total} {stiff} 0.3\n")

        dof = link["DOF"]
        library["link_materials"][name][dof] = mat_total
        mat_total += 1

    for damper in csi.get("LINK PROPERTY DEFINITIONS 04 - DAMPER", []):
        continue
        name = damper["Link"]
        stiff = damper["TransK"]
        dampcoeff = damper["TransC"]
        exp = damper["CExp"]
        model.eval(f"uniaxialMaterial ViscousDamper {mat_total} {stiff} {dampcoeff}' {exp}\n")

        dof = damper["DOF"]
        library["link_materials"][name][dof] = mat_total
        mat_total += 1

    for link in csi.get("LINK PROPERTY DEFINITIONS 10 - PLASTIC (WEN)", []):
        name = link["Link"]

        if not link.get("Nonlinear", False):
            stiff = link["TransKE"]
            model.eval(f"uniaxialMaterial Elastic {mat_total} {stiff}\n")
        else:
            stiff = link["TransK"]
            fy    = link["TransYield"]
            exp   = link["YieldExp"] # TODO
            ratio = link["Ratio"]
            model.eval(f"uniaxialMaterial Steel01 {mat_total} {fy} {stiff} {ratio}\n")

        dof = link["DOF"]
        library["link_materials"][name][dof] = mat_total
        mat_total += 1


    # 2) Frame
    create_frame_sections(csi, model, conv)


    # 3) Shell
    create_shell_sections(csi, model, conv)
    return library


def apply_loads(csi, model):
    "LOAD CASE DEFINITIONS",
    "LOAD PATTERN DEFINITIONS",

    "JOINT LOADS - FORCE",
    "FRAME LOADS - DISTRIBUTED",
    "FRAME LOADS - GRAVITY",
    "FRAME LOADS - POINT",
    "CABLE LOADS - DISTRIBUTED",
    pass



def create_model(csi, types=None, verbose=False):
    """
    Parameters
    ==========
    csi: a dictionary formed by calling ``csi.parse.load("file.b2k")``

    Returns
    =======
    model: opensees.openseespy.Model object
    """

    import opensees.openseespy as ops

    config = CONFIG

    used = {
        "TABLES AUTOMATICALLY SAVED AFTER ANALYSIS"
    }


    #
    # Create model
    #
    dofs = {key:val for key,val in csi["ACTIVE DEGREES OF FREEDOM"][0].items() } # if val }
    dims = {key for key,val in csi["ACTIVE DEGREES OF FREEDOM"][0].items() } # if val }
    ndf = sum(1 for v in csi["ACTIVE DEGREES OF FREEDOM"][0].values())
    ndm = sum(1 for k,v in csi["ACTIVE DEGREES OF FREEDOM"][0].items()
              if k[0] == "U")
    if isinstance(verbose, int) and verbose > 3:
        import sys
        echo_file = sys.stdout
    else:
        echo_file = None
    model = ops.Model(ndm=ndm, ndf=ndf, echo_file=echo_file)

    conv = Converter()

    used.add("ACTIVE DEGREES OF FREEDOM")

#   dofs = [f"U{i}" for i in range(1, ndm+1)]
#   if ndm == 3:
#       dofs = dofs + ["R1", "R2", "R3"]
#   else:
#       dofs = dofs + ["R3"]

    config["ndm"] = ndm
    config["ndf"] = ndf
    config["dofs"] = dofs

    #
    # Create nodes
    #
    create_points(csi, model, None, config, conv)

    # Create materials and sections
    library = create_materials(csi, model, conv)


    # Unimplemented objects
    for item in [
        "CONNECTIVITY - CABLE",
        "CONNECTIVITY - SOLID",
        "CONNECTIVITY - TENDON"]:
        for elem in csi.get(item, []):
            conv.log(UnimplementedInstance(item, elem))

    #
    # Create Links
    #
    create_links(csi, model, library, config, conv)

    #
    # Create frames
    #
    create_frames(csi, model, library, config, conv)

    #
    # Create shells
    #
    create_shells(csi, model, library, conv)

    if verbose and len(conv._log) > 0:
        print_log(conv._log)

    if verbose and False:
        for table in csi:
            if table not in used:
                print(f"\t{table}", file=sys.stderr)

    model.frame_tags = library.get("frame_tags", {})
    return model

