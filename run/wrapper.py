'''
Python class around my MovingHeatSource pybind module.
Might move this properly to the library in the future.
'''
import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
import os
import numpy as np
import meshio
import pdb
import re

class Problem(mhs.Problem):
    # Convenience glob vars
    integerKeys = [
            "nels",
            "maxIter",
            "plot",
            "timeIntegration",
            "sourceTerm",
            "isAdvection",
            ]
    cellMappingMeshio = {
            "line2" : "line",
            "quad4" : "quad",
            "triangle3" : "triangle",
        }
    def __init__(self, caseName):
        self.caseName= caseName
        self.input = {}
        self.iter = 0
        self.postFolder = "post_{}".format( caseName )
        super().__init__()

    # PREPROCESSING
    # Process params
    def readInputFile( self, fileName ):
        self.parseInput( fileName )

    def setMesh( self, points, cells, cell_type ):
        self.input["points"] = points
        self.input["cells"] = cells
        self.input["cell_type"]=cell_type

    def parseInput(self, fileName):
        lines = []

        with open( fileName, 'r') as f:
            lines = f.readlines()

        for idx, l in enumerate(lines):
            lines[idx] =  l.rstrip("\n")

        for l in lines:
            pair = l.split()
            self.input[pair[0]] = float(pair[1])

        for iK in self.integerKeys:
            if iK in self.input:
                self.input[iK] = int(self.input[iK])

    def initialize(self):
        print( "Initializing {}".format( self.caseName ) )
        self.cleanupPrevPost()
        super(Problem, self).initialize( self.input )

    def forceState( self, function ):
        # TODO: Handle prevValues for BDF2+
        val = np.zeros( self.mesh.nnodes )
        for inode in range(self.mesh.nnodes):
            if not(self.mesh.activeNodes[inode]):
                val[inode] = 0
                continue
            pos = self.mesh.posFRF[inode, :]
            val[inode] = function( pos )
        self.initializeIntegrator( val )

    def cleanupPrevPost(self):
        try:
            os.remove( self.caseName + ".pvd" )
            os.removedirs( self.postFolder )
        except OSError:
            pass

    def activate(self, activeElements):
        self.activateDomain( activeElements )

    def iterate(self):
        self.iter += 1
        super(Problem, self).iterate()
        print( "{} iter# {}, time={}".format(
            self.caseName,
            self.iter,
            self.time) )

    def fakeIter(self):
        self.iter += 1
        super(Problem, self).preIterate()

    #POSTPROCESSING
    def writepos( self ):
        os.makedirs(self.postFolder, exist_ok=True)
        #export
        cell_type = self.cellMappingMeshio[self.input["cell_type"]]
        mesh = meshio.Mesh(
            #self.mesh.pos,
            self.mesh.posFRF,
            [ (cell_type, self.mesh.con_CellPoint.con), ],
            point_data={"T": self.unknown.values,
                        "ActiveNodes": self.mesh.activeNodes,},
            cell_data={"ActiveElements":[self.mesh.activeElements]},
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

if __name__=="__main__":
    pass