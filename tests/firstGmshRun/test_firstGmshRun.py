import MovingHeatSource as mhs
from MovingHeatSource.gcode import gcode2laserPath
from CustomStepper import CustomStepper
import numpy as np
import meshio

def deactivateBelowSurface(p,
                           surfacey = 0):
    nels = p.domain.mesh.nels
    activeels = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElement( ielem )
        if (e.getCentroid()[1] < surfacey):
            activeels.append( ielem )
    substrateels = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, activeels )
    p.domain.setActivation( substrateels )

def readMesh(gmshFile):
    m = meshio.read( gmshFile )
    mDict = {}
    mDict["points"] = m.points

    cells = np.vstack( [cell.data for cell in m.cells] )
    mDict["cells"] = cells
    mDict["cell_type"] = "quad4"
    mDict["dimension"] = 2

    return mhs.Mesh( mDict )

def main():
    inputFile = "input.yaml"
    problemInput = mhs.readInput( inputFile )

    # read input
    fixedProblemInput = dict( problemInput )

    # Mesh
    meshFixed = readMesh("untitled.msh")

    pFixed   = mhs.Problem(meshFixed, fixedProblemInput, caseName="fixed")

    # Set path, deactivate
    deactivateBelowSurface( pFixed )

    myDriver = CustomStepper( pFixed,
                              maxAdimtDt=2,
                              threshold=0.025,
                               )

    #while not(pFixed.mhs.path.isOver(pFixed.time)):
    for _ in range(4):
        myDriver.iterate()

def test():
    main()
    reference = "post_fixed_reference/fixed_4.vtu"
    new = "post_fixed/fixed_4.vtu"
    assert mhs.meshio_comparison(reference, new)

if __name__=="__main__":
    test()
