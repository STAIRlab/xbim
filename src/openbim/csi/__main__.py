#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
import math
from openbim.csi import create_model, apply_loads, load, collect_outlines

if __name__ == "__main__":
    import sys

    with open(sys.argv[2], "r") as f:
        csi = load(f)


    if sys.argv[1][1] == "C":
        # Convert
        if "tcl" in sys.argv[1]:
            import xara 
            model = xara.Model(ndm=3, ndf=6, echo_file=sys.stdout)
            model = create_model(csi, model=model, verbose=False)
            sys.exit()
        else:
            model = create_model(csi, verbose=False)
            model.print("-json")
            sys.exit()


    model = create_model(csi, verbose=True)

    if sys.argv[1] == "-E":
        # Eigen
        import veux
        model.constraints("Transformation")
        if len(sys.argv) > 3:
            mode = int(sys.argv[3])
        else:
            mode = 1
        scale = 100
        e = model.eigen(mode)[-1]
        print(f"period = {2*math.pi/math.sqrt(e)}")
        veux.serve(veux.render_mode(model, mode, scale, vertical=3))

    elif sys.argv[1] == "-A":
        # Apply loads and analyze
        apply_loads(csi, model)
        model.analyze(1)

    elif sys.argv[1][:2] == "-V":

        # Visualize
        import veux
        outlines = collect_outlines(csi, model.frame_tags)
        artist = veux.create_artist(model, canvas="gltf", vertical=3,
                    model_config={
                        "frame_outlines": outlines
                    }
        )
        # artist.draw_nodes()
        artist.draw_outlines()
        artist.draw_surfaces()


        if sys.argv[1] == "-Vo":
            artist.save(sys.argv[3])
        else:
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

