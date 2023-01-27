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
            "quad4" : "quad",
            "triangle3" : "triangle",
        }
    def __init__(self, caseName):
        self.caseName= caseName
        self.input = {}
        self.iter = 0
        super().__init__()

    # PREPROCESSING
    # Process params
    def readInputFile( self, fileName ):
        self.parseInput( fileName )

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
        super(Problem, self).initialize( self.input )

    def activate(self, activeElements):
        self.mesh.setActiveElements( activeElements )

    def iterate(self):
        super(Problem, self).iterate()
        self.iter += 1

    #POSTPROCESSING
    def writepos( self, postFolder ):
        os.makedirs(postFolder, exist_ok=True)
        #export
        cell_type = self.cellMappingMeshio[self.input["cell_type"]]
        mesh = meshio.Mesh(
            self.mesh.pos,
            [ (cell_type, self.input["cells"]), ],
            point_data={"T": self.solution},
            cell_data={"ActiveElements":[self.mesh.activeElements]},
        )

        postFilePath = "{}/{}_{}.vtk".format( postFolder, self.caseName, self.iter )
        mesh.write(
            postFilePath,  # str, os.PathLike, or buffer/open file
            # file_format="vtk",  # optional if first argument is a path; inferred from extension
        )

if __name__=="__main__":
    pass
