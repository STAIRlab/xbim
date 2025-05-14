#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
import sys
import veux
import numpy as np
from openbim import inp

if __name__ == "__main__":

    file_name = sys.argv[2]
    lib = inp
    obj = lib.parser.load(file_name, verbose=False)

    model = lib.create_model(obj, verbose=True, mode="visualize")


    if sys.argv[1] == "-C":
        # Convert
        model.print("-json")

    elif sys.argv[1] == "-E":
        # Eigen
        model.constraints("Transformation")
        W = model.eigen(2)
        for w in W:
            print(f"T = {2*np.pi/np.sqrt(w)}")
        veux.serve(veux.render_mode(model, 1, 200.0, vertical=3, canvas="gltf"))

    elif sys.argv[1][:2] == "-A":
        # Apply loads and analyze
        from opensees.helpers import find_nodes, node_average
        top = max(model.nodeCoord(tag)[2] for tag in model.getNodeTags())
        model.pattern("Plain", 1, "Linear")
        for node in find_nodes(model, z=top):
            model.load(node, (0, -1/2, -1, 0, 0, 0))

        # model.integrator("LoadControl", 10.0)
        
        model.algorithm("Linear")
        model.system("UmfPack")
        model.numberer("RCM")
        model.analysis("Static")

        print("Starting analysis")
        model.analyze(1)
        print("Analysis complete")
        if "s" in sys.argv[1]:
            field = {node: stress["sxx"] for node, stress in node_average(model, "stressAtNodes")}
        else:
            field = None
        # artist = veux.render(model, model.nodeDisp, scale=1.0, vertical=3, hide={"node.marker"})
        artist = veux.create_artist(model, vertical=3)
        artist.draw_outlines()
        artist.draw_surfaces(state=model.nodeDisp, field=field, scale=1.0)
        veux.serve(artist)

    elif sys.argv[1] == "-V":
        # Visualize
        import veux
        # veux.serve(veux.render(model, canvas="gltf", vertical=3, hide={"node.marker"}))
                # artist = veux.render(model, model.nodeDisp, scale=1.0, vertical=3, hide={"node.marker"})
        artist = veux.create_artist(model, vertical=3)
        # artist.draw_outlines()
        artist.draw_surfaces()
        veux.serve(artist)

    elif sys.argv[1] == "-Vn":
        # Visualize
        from scipy.linalg import null_space
        model.constraints("Transformation")
        model.analysis("Static")
        K = model.getTangent().T
        v = null_space(K)[:,0] #, rcond=1e-8)
        print(v)


        u = {
            tag: [1000*v[dof-1] for dof in model.nodeDOFs(tag)]
            for tag in model.getNodeTags()
        }

        import veux
        veux.serve(veux.render(model, u, canvas="gltf", vertical=3))

    elif sys.argv[1] == "-Q":
        # Quiet conversion
        pass
    else:
        raise ValueError(f"Unknown operation {sys.argv[1]}")

