#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
import meshio
from . import parser
from .model import create_model


def load(file):
    return meshio.read(file)


def _create_model(obj, verbose=False):
    import opensees.openseespy as ops

    model = ops.Model(ndm=3, ndf=6)

    for i, point in enumerate(obj.points):
        model.node(i, tuple(point))

    # Create Materials

    E = 29e3
    nu = 0.2
    mat = 1
    model.nDMaterial("ElasticIsotropic", mat, E, nu)
    #                                       secTag  E     nu     h    rho
    model.section("ElasticMembranePlateSection", 1, E, 0.25, 1.175, 1.27)

    #
    # Create Elements
    #
    i = 1
    for block in obj.cells:
        if block.type == "hexahedron":
            for nodes in block.data:
                if len(nodes) == 8:
                    model.element("stdBrick", i, tuple(nodes), mat)
                    i += 1
        elif block.type == "quad":
            for nodes in block.data:
                if len(nodes) == 4:
                    model.element("ShellMITC4", i, tuple(nodes), mat, 1)
                    i += 1
                else:
                    print("Quad with ", len(nodes), "nodes")

        elif block.type == "line":
            fsec = 1
            model.section("FrameElastic", fsec, E=E, G=E*0.6, A=1, Iy=1, Iz=1, J=1)
            model.geomTransf("Linear", 1, (0.0, 1.0, 0))
            for nodes in block.data:
                if len(nodes) == 2:
                    model.element("PrismFrame", i, tuple(nodes), section=fsec, transform=1)
                    i += 1

        elif block.type == "triangle":
            for nodes in block.data:
                if len(nodes) == 3:
                    model.element("ShellDKGT", i, tuple(nodes), mat, 1)
                    i += 1
        else:
            print(block.type)


    return model


def apply_loads():
    pass

