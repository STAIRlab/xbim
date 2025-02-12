import warnings 
from .utility import find_row, find_rows
import numpy as np


class _Section:
    def __init__(self, name: str, csi: dict,
                 index: int, model, library, conv):
        self.index = index
        self.name = name
        self.integration = []

        self._create(csi, model, library, conv)

    def _create(self, csi, model, library, conv):
        pass


class _ShellSection(_Section):
    def _create(self, csi, model, library, conv):

        section = find_row(csi["AREA SECTION PROPERTIES"],
                           Section=self.name
        )

        if section is None:
            print(self.name)

        material = find_row(csi["MATERIAL PROPERTIES 01 - GENERAL"],
                            Material=section["Material"]
        )

        material = find_row(csi["MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES"],
                            Material=section["Material"]
        )
        model.section("ElasticMembranePlateSection", self.index,
                      material["E1"],  # E
                      material["G12"]/(2*material["E1"]) - 1, # nu
                      section["Thickness"],
                      material["UnitMass"]
        )
        self.integration.append(self.index)


def create_shell_sections(csi, model, conv):
    tag = 0
    for assign in csi.get("AREA SECTION ASSIGNMENTS", []):
        if assign["Section"] not in library["shell_sections"]:
            library["shell_sections"][assign["Section"]] = \
              _ShellSection(assign["Section"], csi, tag, model, conv)
            tag += len(library["shell_sections"][assign["Section"]].integration)
