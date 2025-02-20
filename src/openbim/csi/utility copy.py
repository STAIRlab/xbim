

class Log:
    def append(self, table, field, object):
        pass

class UnimplementedInstance:
    def __init__(self, name, object, table=None):
        self.name = name 
        self.object = object
        self.table = table

    def __repr__(self):
        return f"{self.name}: {self.object}"

def print_log(log):
    import sys

    types = {i.name for i in log}

    print("Unimplemented features", file=sys.stderr)
    for item in types:
        print(f"\t{item}: {sum(1 for i in log if i.name == item)}", file=sys.stderr)

def find_row(table, **kwds) -> dict:

    for row in table:
        match = True
        for k, v in kwds.items():
            if k not in row or row[k] != v:
                match = False
                break

        if match:
            return row


def find_rows(table, **kwds) -> list:
    rows = []
    for row in table:
        match = True
        for k, v in kwds.items():
            if k not in row or row[k] != v:
                match = False
                break

        if match:
            rows.append(row)

    return rows

