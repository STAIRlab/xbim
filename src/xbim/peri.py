from typing import Dict
import numpy as np

class DimensionMismatch(Exception): pass

class DomainError(Exception): pass

T = np.array

def _midpoint(a: T, b: T, c: T, d: T):
    return 1 / 4 * (a + b + c + d)


def _midpoint(a: T, b: T, c: T, d: T, e: T, f: T, g: T, h: T):
    return 1 / 8 * (a + b + c + d + e + f + g + h)


def _tetvol(a: T, b: T, c: T, d: T):
    v1 = SVector{3}(b - a)
    v2 = SVector{3}(c - a)
    v3 = SVector{3}(d - a)
    vol = 1 / 6 * abs(np.dot(v1, np.cross(v2, v3)))
    return vol


def get_points(nodes: Dict{Int, Vector{Float64}}, elements: Dict{Int, Vector{Int}},
               element_types: Dict{Int, Symbol}):

    nel = len(elements)
    midpoints = np.zeros((3, nel))
    volumes = np.zeros(nel)
    for i in 1:nel
        T = element_types[i]
        if T == Tet4:
            nodeids = elements[i]
            if len(nodeids) == 4:
                a = nodes[nodeids[1]]
                b = nodes[nodeids[2]]
                c = nodes[nodeids[3]]
                d = nodes[nodeids[4]]
                @views midpoints[:, i] = _midpoint(a, b, c, d)
                volumes[i] = _tetvol(a, b, c, d)
            else:
                raise DimensionMismatch("4 nodes needed for element of type Tet4!")


        if T == Hex8:
            nodeids = elements[i]
            if len(nodeids) == 8:
                n1 = nodes[nodeids[1]]
                n2 = nodes[nodeids[2]]
                n3 = nodes[nodeids[3]]
                n4 = nodes[nodeids[4]]
                n5 = nodes[nodeids[5]]
                n6 = nodes[nodeids[6]]
                n7 = nodes[nodeids[7]]
                n8 = nodes[nodeids[8]]
                midpoints[:, i] = _midpoint(n1, n2, n3, n4, n5, n6, n7, n8)
                volume1 = _tetvol(n1, n2, n4, n5)
                volume2 = _tetvol(n2, n3, n4, n7)
                volume3 = _tetvol(n2, n5, n6, n7)
                volume4 = _tetvol(n4, n5, n7, n8)
                volume5 = _tetvol(n2, n4, n5, n7)
                volumes[i] = volume1 + volume2 + volume3 + volume4 + volume5
            else:
                raise DimensionMismatch("8 nodes needed for element of type Hex8")

        else:
            msg = "Element of type $T not supported!\n"
            raise DomainError(T, msg)


    return midpoints, volumes

