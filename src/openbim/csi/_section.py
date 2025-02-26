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


def _create_shell_integration(csi, model, conv):
    pass

def _create_shell_section(csi, assign, model, conv):

    section = find_row(csi["AREA SECTION PROPERTIES"],
                        Section=assign["Section"]
    )

    if section is None:
        # TODO: log
        print(assign["Section"])

    tag = conv.define("ShellSection", "section", assign["Section"])

    material = find_row(csi["MATERIAL PROPERTIES 01 - GENERAL"],
                        Material=section["Material"]
    )

    material = find_row(csi["MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES"],
                        Material=section["Material"]
    )
    model.section("ElasticShell", tag,
                    material["E1"],  # E
                    material["G12"]/(2*material["E1"]) - 1, # nu
                    section["Thickness"],
                    material["UnitMass"]
    )
    # self.integration.append(self.index)
    return tag


def create_shell_sections(csi, model, conv):
    tag = 0
    for assign in csi.get("AREA SECTION ASSIGNMENTS", []):
        if not conv.identify("ShellSection", "section",  assign["Section"]):

            if not _create_shell_section(csi, assign, model, conv):
                warnings.warn(f"Section {assign['Section']} not found")
                continue
            # tag += len(library["shell_sections"][assign["Section"]].integration)

