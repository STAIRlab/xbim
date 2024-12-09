#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
import json
import shlex

CONSTANTS = {
        "Yes": True,
        "No":  False
}

def _parse_value(v):
    if v in CONSTANTS:
        return CONSTANTS[v]

    try:
        return json.loads(v)
    except:
        return v

def load(file, append: dict=None):
    """
    Read file-like object file and form a dictionary.
    """
    if append is None:
        tables = {}
    else:
        tables = append

    current_table = None
    current_item  = None
    for line in file:
        if "END TABLE DATA" in line:
            break

        # Skip empty lines
        if line.isspace():
            continue

        if "TABLE:" in line:
            table_name = shlex.split(line)[1]
            current_item  = None

            # Append if table exists (append argument given)
            if table_name in tables:
                current_table = tables[table_name]
            else:
                current_table = []

            tables[table_name] = current_table


        # Data line
        elif current_table is not None:

            if current_item is None:

                current_item = {}
                current_table.append(current_item)

            continue_item = False
            tokens = shlex.split(line)
            for i,kv in enumerate(tokens):
                if kv == "_":
                    if i == len(tokens)-1:
                        continue_item = True
                        break

                    # Sometimes there is a random "_" in the middle
                    # of a line?
                    else:
                        continue

                k, v = kv.split("=", maxsplit=1)
                current_item[k] = _parse_value(v)

            if not continue_item:
                current_item = None

    return tables

