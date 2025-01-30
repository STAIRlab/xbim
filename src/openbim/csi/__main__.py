from openbim.csi import create_model, apply_loads, load, collect_outlines

if __name__ == "__main__":
    import sys

    with open(sys.argv[2], "r") as f:
        csi = load(f)


    model = create_model(csi, verbose=True)

    if sys.argv[1] == "-C":
        # Convert
        model.print("-json")

    elif sys.argv[1] == "-E":
        # Eigen
        import veux
        model.constraints("Transformation")
        model.eigen(2)
        veux.serve(veux.render_mode(model, 1, 50.0, vertical=3))

    elif sys.argv[1] == "-A":
        # Apply loads and analyze
        apply_loads(csi, model)
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

