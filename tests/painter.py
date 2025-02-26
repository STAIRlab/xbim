import numpy as np
import os
import io
from openbim.csi import create_model, apply_loads, load, collect_outlines

if __name__ == "__main__":
    import sys

    model_file = "models/csi/painter/Painter_Street_v1.2.b2k"
    with open(model_file, "r") as f:
        csi = load(f)

    # This is an opensees.openseespy.Model class
    # https://opensees.stairlab.io/user/manual/model/model_class.html
    model = create_model(csi, verbose=True)

    if sys.argv[1] == "-C":
        # Convert
        model.print("-json")

    elif sys.argv[1] == "-E":
        # Eigen analysis
        import veux
        if len(sys.argv) > 2:
            mode = int(sys.argv[2])
        else:
            mode = 1
        scale = 25

        #
        model.constraints("Transformation")
        e = model.eigen(mode)[-1]
        print(f"period = {2*np.pi/np.sqrt(e)}")

        # documentation on the canvas argument is here: https://veux.io/library/canvas.html
        artist = veux.render_mode(model, mode, scale, vertical=3, canvas="plotly")
        veux.serve(artist)

    elif sys.argv[1] == "-A":
        # Apply loads and analyze
        apply_loads(csi, model)
        model.analysis("Static")
        model.analyze(1)

    elif sys.argv[1][:2] == "-V":
        # Visualize

        import veux
        outlines = collect_outlines(csi, model.frame_tags)

        artist = veux.render(model, canvas="gltf", vertical=3,
                    reference={"frame.surface", "frame.axes"},
                    model_config={
                        "frame_outlines": outlines
                    }
        )

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

