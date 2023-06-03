'''
Python class around my MovingHeatSource pybind module.
Might move this properly to the library in the future.
'''
import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Release/')
import MovingHeatSource as mhs
import os, shutil
import numpy as np
import meshio
import pdb
import re

def readInput(fileName):
    integerKeys = [
            "nels",
            "maxIter",
            "plot",
            "timeIntegration",
            "sourceTerm",
            "isAdvection",
            ]
    problemInput = {}
    lines = []
    with open( fileName, 'r') as f:
        lines = f.readlines()

    for idx, l in enumerate(lines):
        lines[idx] =  l.rstrip("\n")

    for l in lines:
        pair = l.split()
        problemInput[pair[0]] = float(pair[1])

    for iK in integerKeys:
        if iK in problemInput:
            problemInput[iK] = int(problemInput[iK])

    return problemInput

class Problem(mhs.Problem):
    # Convenience glob vars
    cellMappingMeshio = {
            "line2" : "line",
            "quad4" : "quad",
            "triangle3" : "triangle",
        }
    def __init__(self, mesh, input, problem=None, caseName="case"):
        self.caseName= caseName
        self.input = input
        self.iter = 0
        self.postFolder = "post_{}".format( self.caseName )

        if problem:
            self.input = problem.input
            super().__init__(problem)
            self.setPointers()
            self.iter = problem.iter
        else:
            super().__init__(mesh, input)

    # PREPROCESSING
    # Process params
    def setMesh( self, points, cells, cell_type ):
        self.input["points"] = points
        self.input["cells"] = cells
        self.input["cell_type"]=cell_type

    def forceState( self, function ):
        # TODO: Handle prevValues for BDF2+
        val = np.zeros( self.domain.mesh.nnodes )
        for inode in range(self.domain.mesh.nnodes):
            if not(self.domain.activeNodes.x[inode]):
                val[inode] = 0
                continue
            pos = self.domain.mesh.posFRF[inode, :]
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

    def preiterate(self, canPreassemble):
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

    #POSTPROCESSING
    def writepos( self,
                 rf = "FRF",
                 shift=None,
                 functions={},
                 nodeMeshTags={},
                 cellMeshTags={},
                 ):
        '''
        rf : reference frame, either FRF or MRF
        '''
        os.makedirs(self.postFolder, exist_ok=True)

        if   rf=="FRF":
            pos = self.domain.mesh.posFRF
        elif rf=="MRF":
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
        cellType = self.cellMappingMeshio[myCellType]
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
        timestepLine = '<DataSet timestep="{}" group="" part="0" file="{}"/>\n'.format(round(self.time, 3), postFilePath)

        pvdContents = []
        with open(pvdFileName, "r") as f:
            pvdContents = f.readlines()

        pvdContents.insert(-2, timestepLine)

        with open(pvdFileName, "w") as f:
            pvdContents = "".join(pvdContents)
            f.write(pvdContents)


def meshio_comparison(ref, new, tol=1e-7):
    # Load datasets
    refds = meshio.read( ref )
    newds = meshio.read( new )
    # Compare mesh
    # Compare points
    pdiff = np.abs( refds.points - newds.points )
    if not( (pdiff < tol).all() ):
        return False
    # Compare connectivities
    for ref_cblock, new_cblock in zip( refds.cells, newds.cells ):
        if not( (ref_cblock.data == new_cblock.data).all() ):
            return False
    # Compare point data
    for key in refds.point_data.keys():
        refpdata = refds.point_data[key]
        newpdata = newds.point_data[key]
        pdatadiff = np.abs(refpdata - newpdata )
        if not( (pdatadiff < tol).all() ):
            return False
    # Compare cell data
    for key in refds.cell_data.keys():
        for refcdata, newcdata in zip( refds.cell_data[key], newds.cell_data[key] ):
            cdatadiff = np.abs( refcdata - newcdata )
            if not( (cdatadiff < tol).all() ):
                return False
    # All comparisons passed
    return True

if __name__=="__main__":
    pass
