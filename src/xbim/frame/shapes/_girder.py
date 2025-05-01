
class GirderSection(Shape):
    def __init__(self,
        thickness_top  : float,
        thickness_bot  : float,
        height         : float,
        width_top      : float,
        width_webs     : list,
        web_spacing    : float,
        web_slope      : float = 0.0,
        overhang       : float = None
        ):
        self._thickness_top  = thickness_top
        self._thickness_bot  = thickness_bot
        self._height         = height
        self._width_top      = width_top
        self._width_webs     = width_webs
        self._web_spacing    = web_spacing
        self._web_slope      = web_slope
        self._overhang       = overhang

    def create_mesh(self):
        #                                ^ y
        #                                |
        # _  |_______________________________________________________|
        #    |_____  _______________ _________ _______________  _____|
        #          \ \             | |       | |             / /
        #           \ \            | |   |   | |            / /
        #            \ \___________| |_______| |___________/ /
        # _           \__________________+__________________/  ---> x
        #             |                                     |

        import opensees.units
        spacing = opensees.units.units.spacing
        thickness_top  = self._thickness_top
        thickness_bot  = self._thickness_bot
        height         = self._height
        width_top      = self._width_top
        width_webs     = self._width_webs
        web_spacing    = self._web_spacing
        web_slope      = self._web_slope
        overhang       = self._overhang

        # Dimensions
        #------------------------------------------------------
        inside_height = height - thickness_bot - thickness_top


        # width of bottom flange
        if overhang:
            width_bot = width_top - \
                    2*(overhang + web_slope*(inside_height + thickness_bot))
        else:
            width_bot = web_centers[-1] - web_centers[0] \
                    + width_webs[1]/2 + width_webs[0]/2

        # number of internal web *spaces*
        niws = len(width_webs) - 3

        # list of web centerlines?
        web_centers   = [
            -width_bot/2 - inside_height/2*web_slope + 0.5*width_webs[1],
            *niws @ spacing(web_spacing, "centered"),
            width_bot/2 + inside_height/2*web_slope - 0.5*width_webs[-1]
        ]

        # Build section
        #------------------------------------------------------
        girder_section = [
            # add rectangle patch for top flange
            patch.rect(corners=[
                [-width_top/2, height - thickness_top],
                [+width_top/2, height                ]]),

            # add rectangle patch for bottom flange
            patch.rect(corners=[
                [-width_bot/2,        0.0      ],
                [+width_bot/2,  +thickness_bot]]),

            # sloped outer webs
            patch.rhom(
                height = inside_height,
                width  = width_webs[0],
                slope  = -web_slope,
                center = [web_centers[0], thickness_bot + inside_height/2]
            ),
            patch.rhom(
                height = inside_height,
                width  = width_webs[-1],
                slope  = web_slope,
                center = [web_centers[-1], thickness_bot + inside_height/2]
            )
        ] + [
            patch.rect(corners=[
                [loc - width/2,        thickness_bot],
                [loc + width/2,  height - thickness_top]]
            )
            for width, loc in zip(width_webs[1:-1], web_centers[1:-1])
        ]

        return create_mesh(girder_section, mesh_size=min(thickness_bot, thickness_top, *width_webs)/3.0)


