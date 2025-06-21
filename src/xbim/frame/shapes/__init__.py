from xsection.library import (
    WideFlange,
    Rectangle,
    HollowRectangle,
    Angle,
    Channel,
    HollowRectangle,
    from_aisc
)

if __name__ == "__main__":
    import veux
    d  = 100
    tw = 3
    bf = 75
    tf = 3

    mesh = WideFlange(d=d, b=bf, t=tf, tw=tw).create_shape()

    print(mesh.summary())

#   from shps.frame.solvers.plastic import PlasticLocus
#   PlasticLocus(mesh).plot()#(phi=0.5, ip=5)
#   import matplotlib.pyplot as plt
#   plt.show()

    artist = veux.create_artist(((mesh.mesh.nodes, mesh.mesh.cells())), ndf=1)

    field = mesh.torsion.warping()
    field = {node: value for node, value in enumerate(field)}

    artist.draw_surfaces(field = field)
    artist.draw_origin()
#   R = artist._plot_rotation
#   artist.canvas.plot_vectors([R@[*geometry.centroid, 0] for i in range(3)], R.T)
    artist.draw_outlines()
    veux.serve(artist)
