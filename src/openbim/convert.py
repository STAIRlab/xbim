import re
from collections import defaultdict

RE = {
    "joint_key": re.compile("Joint[0-9]")
}

TYPES = {
    "Shell": {
        "Elastic": "ShellMITC4",
    },
    "Frame": {
        "Elastic": "PrismFrame"
    }
}

class Converter:
    def __init__(self):
        self._log = []
        self._csi_names = {
            # Abaqus
            "Material": {
                "material": {},
                "uniaxial": {}
            },
            # CSI
            "Joint": {
                "node": {
                    # 1:    1,
                    # "P2": 2,
                    # 3:    3
                },
            },
            "Shell": {"element": {}},
            "Link":  {"element": {}},
            "Frame": {"element": {}},
            "AnalSect": {
                "section": {},
                "integration": {}
            },
            "ShellSection": {
                "section": {}
            }
        }
        self._ops_count = {
                "node":       0,
                "transform":  0,
                "element":    0,
                "nDMaterial": 0,
                "section":    0,
                "material":   0,
                "integration": 0
        }

        self._library = {
            "frame_sections": {},
            "shell_sections": {},
            "link_materials": defaultdict(dict),
        }

    def identify(self, csi_type, ops_type, csi_name)->int:
        if ops_type not in self._csi_names[csi_type]:
            self._csi_names[csi_type][ops_type] = {}
            return None

        if csi_name in self._csi_names[csi_type][ops_type]:
            return self._csi_names[csi_type][ops_type][csi_name]

        return None

    def define(self, csi_type, ops_type, csi_name=None, item=None)->int:
        if ops_type not in self._csi_names[csi_type]:
            self._csi_names[csi_type][ops_type] = {}

        if csi_name is None:
            self._ops_count[ops_type] += 1
            return self._ops_count[ops_type]

        if csi_name not in self._csi_names[csi_type][ops_type]:
            self._ops_count[ops_type] += 1
            if item is None:
                item = self._ops_count[ops_type]
            self._csi_names[csi_type][ops_type][csi_name] = item

        return self._csi_names[csi_type][ops_type][csi_name]


    def log(self, message):
        self._log.append(message)

