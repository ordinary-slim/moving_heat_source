'''
Python class around my MovingHeatSource pybind module.
'''
from MovingHeatSource.cpp import *
import os, shutil
import numpy as np
import meshio
import yaml
from MovingHeatSource.gcode import gcode2laserPath

def readInput(fileName):
    params = {}
    with open(fileName, 'r') as paramsFile:
        params = yaml.safe_load(paramsFile)
    for k, v in params.items():
        if isinstance( v, list ):
            params[k] = np.array( v )
    return params

class Problem(Problem):
    # Convenience glob vars
    def __init__(self, mesh, input, problem=None, caseName="case"):
        self.caseName= caseName
        self.input = dict(input)
        self.iter = 0
        self.postFolder = "post_{}".format( self.caseName )

        if problem:
            self.input = problem.input
            super().__init__(problem)
            self.setPointers()
            self.iter = problem.iter
        else:
            super().__init__(mesh, input)
        self.cleanupPrevPost()

    # PREPROCESSING
    def forceState( self, function ):
        # TODO: Handle prevValues for BDF2+
        val = np.zeros( self.domain.mesh.nnodes )
        for inode in range(self.domain.mesh.nnodes):
            if not(self.domain.activeNodes.x[inode]):
                val[inode] = 0
                continue
            pos = self.domain.posLab[inode, :]
            val[inode] = function( pos )
        self.initializeIntegrator( val )

    def cleanupPrevPost(self):
        try:
            os.remove( self.caseName + ".pvd" )
        except FileNotFoundError:
            pass

        try:
            shutil.rmtree( self.postFolder )
        except FileNotFoundError:
            pass

    def preiterate(self, canPreassemble=True):
        self.iter += 1
        super(Problem, self).preIterate(canPreassemble)

    def iterate(self):
        if not(self.hasPreIterated):
            self.preiterate(True)
        super(Problem, self).iterate()
        print( "{} iter# {}, time={}".format(
            self.caseName,
            self.iter,
            self.time) )

    def fakeIter(self):
        self.iter += 1
        super(Problem, self).preIterate(True)
        print( "{} fake iter# {}, time={}".format(
            self.caseName,
            self.iter,
            self.time) )

    def frf2mrf(self, speed=None):
        if speed is None:
            speed = self.mhs.speed
        self.domain.mesh.setSpeedFRF( speed )
        self.setAdvectionSpeed( -speed )
        self.mhs.setSpeed( np.zeros(3) )

    def setPath(self, gcodeFile):
        self.mhs.setPath( *gcode2laserPath( gcodeFile, defaultPower=self.input["power"] ) )

    #POSTPROCESSING
    def writepos( self,
                 rf = "Lab",
                 shift=None,
                 functions={},
                 nodeMeshTags={},
                 cellMeshTags={},
                 ):
        '''
        rf : reference frame, either HeatSource or Lab
        '''
        os.makedirs(self.postFolder, exist_ok=True)

        if   rf=="Lab":
            pos = self.domain.posLab
        elif rf=="HeatSource":
            pos = self.domain.mesh.pos
        else:
            print("Wrong value of rf")
            exit()

        # Debugging feature for correctness checks
        if shift is not None:
            pos = np.array( pos )
            for irow in range(pos.shape[0]):
                pos[irow,:] += shift

        point_data={"T": self.unknown.values,
                    "Pulse": self.mhs.pulse,
                    "ActiveNodes": self.domain.activeNodes.x,
                    }
        cell_data={"ActiveElements":[self.domain.activeElements.x]}
        for label, fun in functions.items():
            point_data[label] = fun.values

        for label, tag in nodeMeshTags.items():
            point_data[label] = tag.x
        for label, tag in cellMeshTags.items():
            cell_data[label] = [tag.x]

        #export
        myCellType = self.domain.mesh.elementTypes[0].name
        cellMappingMeshio = {
                "line2" : "line",
                "quad4" : "quad",
                "triangle3" : "triangle",
                "hexa8" : "hexahedron",
            }
        cellType = cellMappingMeshio[myCellType]
        mesh = meshio.Mesh(
            pos,
            [ (cellType, self.domain.mesh.con_CellPoint.con), ],
            point_data=point_data,
            cell_data=cell_data,
        )

        postFilePath = "{}/{}_{}.vtu".format( self.postFolder, self.caseName, self.iter )
        mesh.write(
            postFilePath,  # str, os.PathLike, or buffer/open file
        )
        #update pvd
        pvdFileName =  self.caseName + ".pvd" 
        if not( os.path.isfile(pvdFileName) ):
            baseLinesPVD = [
                '<?xml version="1.0"?>',
                '<VTKFile type="Collection" version="0.1">',
                '<Collection>',
                '</Collection>',
                '</VTKFile>',
                ]
            with open( pvdFileName, 'w' ) as pvd:
                pvd.write( "\n".join(baseLinesPVD) )
        #insert time-step into pvd.
        #always two lines before last
        timestepLine = '<DataSet timestep="{:.6f}" group="" part="0" file="{}"/>\n'.format(self.time, postFilePath)

        pvdContents = []
        with open(pvdFileName, "r") as f:
            pvdContents = f.readlines()

        pvdContents.insert(-2, timestepLine)

        with open(pvdFileName, "w") as f:
            pvdContents = "".join(pvdContents)
            f.write(pvdContents)


def meshio_comparison(ref, new,
                      psets=None,
                      csets=None,
                      tol=1e-7):
    isSame = True
    differentKeys = []
    # Load datasets
    refds = meshio.read( ref )
    newds = meshio.read( new )
    # Compare mesh
    # Compare points
    pdiff = np.abs( refds.points - newds.points )
    if not( (pdiff < tol).all() ):
        isSame = False
        differentKeys.append("Points")
    # Compare connectivities
    for ref_cblock, new_cblock in zip( refds.cells, newds.cells ):
        if not( (ref_cblock.data == new_cblock.data).all() ):
            isSame = False
            differentKeys.append("Connectivity")
    # Compare point data
    if not psets:
        try:
            len(psets)
        except TypeError:
            psets = refds.point_data.keys()
    for key in psets:
        keyIsSame = True
        refpdata = refds.point_data[key]
        newpdata = newds.point_data[key]
        pdatadiff = np.abs(refpdata - newpdata )
        if not( (pdatadiff < tol).all() ):
            keyIsSame = False
        if not(keyIsSame):
            isSame = False
            differentKeys.append( key )
    # Compare cell data
    if not csets:
        try:
            len(csets)
        except TypeError:
            csets = refds.cell_data.keys()
    for key in csets:
        keyIsSame = True
        for refcdata, newcdata in zip( refds.cell_data[key], newds.cell_data[key] ):
            cdatadiff = np.abs( refcdata - newcdata )
            if not( (cdatadiff < tol).all() ):
                keyIsSame = False
            if not(keyIsSame):
                isSame = False
                differentKeys.append( key )
    # All comparisons passed
    if not(isSame):
        print("{}:\nTHESE KEYS ARE NOT MATCHING: {}".format( new, differentKeys) )
    return isSame

try:
    import gmsh
    _has_gmsh = True
except ModuleNotFoundError:
    _has_gmsh = False

if _has_gmsh:
    def gmshModelToMesh( model: gmsh.model ):
        '''
        gmsh model to points and cells
        Taken from Dolfinx project:
        https://github.com/FEniCS/dolfinx/blob/main/python/dolfinx/io/gmshio.py
        Doesn't support mixed meshes nor multiple physical groups
        '''
        gmshElType2myElType = {
                1 : "line2",
                2 : "triangle3",
                3 : "quad4",
                5 : "hexa8",
                }
        dim = model.getDimension()
        # POINTS & INDICES
        indices, points, _ = model.mesh.getNodes()
        points = points.reshape(-1, 3)

        #Gmsh indices starts at 1
        indices -= 1

        # In some cases, Gmsh does not return the points in the same
        # order as their unique node index.
        # We therefore sort nodes in
        # geometry according to the unique index
        perm_sort = np.argsort(indices)
        assert np.all(indices[perm_sort] == np.arange(len(indices)))
        points = points[perm_sort]
        if dim < 3:
            keepDim = [True]*3
            for idx in range(3):
                if (np.max(points[:,idx]) - np.min(points[:,idx])) < 1e-7:
                    keepDim[idx] = False
            points = np.transpose( np.vstack( tuple( [points[:, idx] for idx in range(3) if keepDim[idx]] ) ) )


        # CONNECTIVITY
        domainCellType = None
        domainConnectivity = None
        # Only extract connectivity from max dim of model
        _, domainTag = model.getPhysicalGroups(dim=dim)[0]#only first physical group
        subdomains = model.getEntitiesForPhysicalGroup( dim, domainTag )
        subdomainConnectivities = []
        for domain in subdomains:
            (entity_types, entity_tags, entity_connectivities) = model.mesh.getElements(dim, tag=domain)
            assert len(entity_types) == 1# no support for multiple cell types
            # Determine number of local nodes per element to create the
            # topology of the elements
            properties = model.mesh.getElementProperties(entity_types[0])
            name, dim, _, num_nodes, _, _ = properties

            subdomainConnectivity = entity_connectivities[0].reshape(-1, num_nodes) - 1
            subdomainConnectivities.append( (entity_types[0], subdomainConnectivity) )
        # Condense all connectivities into one
        domainCellType = subdomainConnectivities[0][0]
        domainConnectivity = subdomainConnectivities[0][1]
        for subdomainCellType, subdomainConnecitivity in subdomainConnectivities[1:]:
            if not(subdomainCellType == domainCellType):
                raise ValueError("Mixed element meshes are not supported")
            domainConnectivity = np.concatenate( (domainConnectivity, subdomainConnecitivity), axis=0 )

        return Mesh( {"points":points,
                      "cells":domainConnectivity,
                      "cell_type":gmshElType2myElType[ domainCellType ],
                      "dimension":model.getDimension(), }
                    )
    

if __name__=="__main__":
    pass
