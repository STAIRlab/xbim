import sys
import ifcopenshell

import openbim.msh as g2o
import opensees.openseespy as ops 

import  export 
import meshing
import analysis



ifcFile = ifcopenshell.open(sys.argv[1])
elements = ifcFile.by_type('IfcElement')


data = export. exportProperties(ifcFile, elements)

FileName = 'Facade'
file = FileName + '.step'


print("STEPwriter")
stepfile = export.STEPwriter(elements, data, FileName)

print("mesh_physical_groups")
MyModel = meshing.mesh_physical_groups(stepfile, data, True)


print("fix_boundaries")
Fixed = meshing.fix_boundaries(MyModel)

print("meshing")
Meshed = meshing.meshing(Fixed, True)

print("Create tets")
#model, elementTags, finalTags, nodeTags
opsModel = analysis.Create4NodesTetraedron(Meshed, data)

print("Static analysis")
staticAnalysis = analysis.StaticAnalysis(opsModel[0], opsModel[2], Meshed )

EigenAnalysis = analysis.EigenValue(opsModel[0], Meshed )

