from xsection.library import (
    WideFlange,
    Rectangle,
    HollowRectangle,
    Angle,
    Channel,
    HollowRectangle,
    from_aisc
)

if False:
    def from_aisc(identifier: str, mesh:tuple=None, units=None, ndm=None, tag=None, **kwds):
        if units is None:
            import opensees.units.english as _units
            units = _units

        aisc_data = load_aisc(identifier, units=units)

        if aisc_data is None:
            raise ValueError(f"Cannot find section with identifier {identifier}")


        if identifier[0] == "W":
            return WideFlange(d=aisc_data["d"],
                            b=aisc_data["bf"],
                            t=aisc_data["tf"],
                            tw=aisc_data["tw"],
                            k=aisc_data.get("k1", None),
                            **kwds)

        elif identifier[0] == "L":
            return Angle(d=aisc_data["d"],
                        b=aisc_data["b2"],
                        t=aisc_data["t"])

        elif identifier[0] == "C":
            return Channel(d=aisc_data["d"],
                        b=aisc_data["bf"],
                        tf=aisc_data["tf"],
                        tw=aisc_data["tw"])
        elif identifier[:3] == "HSS":
            return HollowRectangle(d=aisc_data["h/tdes"]*aisc_data["tdes"],
                                b=aisc_data["b/tdes"]*aisc_data["tdes"],
                                t=aisc_data["tdes"],
            )
        else:
            raise ValueError(f"Cannot create section from identifier {identifier}")



    def load_aisc(SectionName, props="", units=None)->dict:
        """Load cross section properties from AISC database.

        props:
            A list of AISC properties, or one of the following:
            - 'simple': `A`, `Ix`, `Zx`

        """

        if units is None:
            import opensees.units.english as _units
            units = _units

        from shps.frame.shapes.aisc_imperial import imperial
        SectData = imperial[SectionName.upper()]


        if props == "simple":
            props = ""
            return

        elif props:
            props = props.replace(" ", "").split(",")
            sectData = {k: v for k, v in SectData.items() if k in props}
            if "I" in props:
                sectData.update({"I": SectData["Ix"]})
            return sectData

        for k,v in list(SectData.items()):
            try:
                SectData[k] = float(v)
            except:
                continue


        UNITS = [
            ("d"  ,  units.inch ),
            ("k1"  , units.inch ),
            ("h/tdes", 1 ),
            ("b/tdes", 1 ),
            ("tdes", units.inch ),
            ("bf" , units.inch ),
            ("tw" , units.inch ),
            ("tf" , units.inch ),
            ("A"  , units.inch**2 ),
            ("Ix" , units.inch**4 ),
            ("Iy" , units.inch**4 ),
        ]

        return {
            k: SectData[k]*scale
            for k, scale in UNITS if k in SectData and isinstance(SectData[k], (float, int))
        }



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
