import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Release/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
from wrapper import Problem
import pdb

def mesh(box):
    cell_type="quad4"
    meshDen = 4
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
        np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( int )
    return points, cells, cell_type

def isInsideBox( mesh, box ):
    activeElements = []
    for ielem in range( mesh.nels ):
        el = mesh.getElement( ielem )
        pos = mesh.posFRF[el.con]
        xmin = min(pos[:, 0])
        xmax = max(pos[:, 0])
        ymin = min(pos[:, 1])
        ymax = max(pos[:, 1])

        isInside = 1*(xmin>=box[0] and xmax <= box[1] and ymin >= box[2] and ymax <= box[3])
        activeElements.append(isInside)

    activeElements = mhs.MeshTag( mesh, mesh.dim, activeElements )
    return activeElements

def setAdimR( adimR, p ):
    r = p.input["radius"]
    HeatSourceSpeedX = max( abs(p.input["HeatSourceSpeedX"]), abs(p.input["advectionSpeedX"]))
    HeatSourceSpeedY = max( abs(p.input["HeatSourceSpeedY"]), abs(p.input["advectionSpeedY"]))
    HeatSourceSpeedZ = max( abs(p.input["HeatSourceSpeedZ"]), abs(p.input["advectionSpeedZ"]))
    speed  = np.linalg.norm( np.array( [HeatSourceSpeedX, HeatSourceSpeedY, HeatSourceSpeedZ] ) )
    return (adimR * r / speed)

def getMaxT( p ):
    maxT = -1
    posMaxT = None
    for inode in range(p.domain.mesh.nnodes):
        if (p.unknown.values[inode] > maxT):
            maxT = p.unknown.values[inode]
            posMaxT = p.domain.mesh.pos[inode]
    return maxT, posMaxT

def debugHeatSourceNPeak( p ):
    print("Position on heat source in Xi:", p.mhs.currentPosition)
    print("MaxT = {}, pos max T Xi = {}".format( *getMaxT( p ) ))

if __name__=="__main__":
    inputFile = "input.txt"
    boxRef = [-16, 16, -5, 5]
    boxInac = [-32, 32, -5, 5]
    adimR = 1

    pFineFRF         = Problem("fineFRF")
    pFRF             = Problem("FRF")
    pNoTransportMRF           = Problem("NoTransportMRF")
    pTransportedMRF             = Problem("TransportedMRF")
    pMRFTransporter  = Problem("MRFTransporter")

    meshFineFRF       = mhs.Mesh()
    meshFRF           = mhs.Mesh()
    meshNoTransportMRF= mhs.Mesh()
    meshTransportedMRF= mhs.Mesh()
    meshMRFTransporter= mhs.Mesh()



    meshPhys = mhs.Mesh()
    meshInputPhys = {}
    meshInputPhys["points"], meshInputPhys["cells"], meshInputPhys["cell_type"] = mesh(boxRef)

    meshBg = mhs.Mesh()
    meshInputBg = {}
    meshInputBg["points"], meshInputBg["cells"], meshInputBg["cell_type"] = mesh(boxInac)

    #read input
    for p in [pFineFRF, pFRF, pNoTransportMRF, pTransportedMRF, pMRFTransporter]:
        p.parseInput( inputFile )
        p.input["cell_type"] = meshInputPhys["cell_type"]#TODO: quickfix!

    # set dt
    dt = setAdimR( adimR, pFRF )
    for p in [pFRF, pNoTransportMRF, pTransportedMRF, pMRFTransporter]:
        p.input["dt"] = dt
    ##determine fine problem tstep size
    approxFine_dt = pow(dt, 2)
    approxFine_dt = min( approxFine_dt, dt / 32.0 )
    fineStepsPerStep = int( np.ceil( dt / approxFine_dt ) )
    fine_dt = dt / float( fineStepsPerStep )
    pFineFRF.input["dt"] = fine_dt

    #set MRF business NO TRANSPORT
    for p in [pNoTransportMRF,]:
        p.input["isAdvection"] = 1
        p.input["advectionSpeedX"] = -pTransportedMRF.input["HeatSourceSpeedX"]
        p.input["speedFRF_X"]      = pTransportedMRF.input["HeatSourceSpeedX"]
        p.input["HeatSourceSpeedX"] = 0.0
    #set MRF business TRANSPORT
    for p in [pTransportedMRF]:
        p.input["isAdvection"] = 1
        p.input["advectionSpeedX"] = -pTransportedMRF.input["HeatSourceSpeedX"]
        p.input["speedFRF_X"]      = pTransportedMRF.input["HeatSourceSpeedX"]
        p.input["HeatSourceSpeedX"] = 0.0

    for p, m in zip([pFineFRF, pFRF, pMRFTransporter], [meshFineFRF, meshFRF, meshMRFTransporter]):
        m.initializeMesh( meshInputPhys )
        p.initialize(m)
    for p, m in zip([pNoTransportMRF, pTransportedMRF,], [meshNoTransportMRF, meshTransportedMRF,]):
        m.initializeMesh( meshInputBg )
        p.initialize(m)

    pMRFTransporter.unknown.interpolate(  pTransportedMRF.unknown )
    for tF, sF in zip(pMRFTransporter.previousValues, pTransportedMRF.previousValues):
        tF.interpolate( sF )

    maxIter = pFRF.input["maxIter"]
    # FORWARD
    for iteration in range(maxIter):
        #fine problem
        '''
        for istep in range(fineStepsPerStep):
            pFineFRF.iterate()
        pFineFRF.writepos()
        '''

        #iter FRF
        pFRF.iterate()#assembly + solve
        pFRF.writepos()

        #iter NoTransportMRF
        pNoTransportMRF.updateFRFpos()
        activeElements = isInsideBox( pNoTransportMRF.domain.mesh, boxRef )
        pNoTransportMRF.domain.setActivation( activeElements )
        print("---BEFORE-----------")
        debugHeatSourceNPeak( pNoTransportMRF )
        print("--------------------")

        pNoTransportMRF.iterate()
        print("---AFTER------------")
        debugHeatSourceNPeak( pNoTransportMRF )
        print("--------------------")
        #pdb.set_trace()

        pNoTransportMRF.writepos()

        #iter transportedMRF
        pTransportedMRF.updateFRFpos()
        pTransportedMRF.unknown.interpolate( pMRFTransporter.unknown )
        for tF, sF in zip(pTransportedMRF.previousValues, pMRFTransporter.previousValues):
            tF.interpolate( sF )

        activeElements = isInsideBox( pTransportedMRF.domain.mesh, boxRef )
        pTransportedMRF.domain.setActivation( activeElements )
        print("---BEFORE-----------")
        debugHeatSourceNPeak( pTransportedMRF )
        print("--------------------")

        pTransportedMRF.iterate()
        print("---AFTER------------")
        debugHeatSourceNPeak( pTransportedMRF )
        print("--------------------")
        #pdb.set_trace()

        pMRFTransporter.fakeIter()
        pMRFTransporter.unknown.interpolate( pTransportedMRF.unknown )
        for tF, sF in zip(pMRFTransporter.previousValues, pTransportedMRF.previousValues):
            tF.interpolate( sF )
        pTransportedMRF.writepos()
        pMRFTransporter.writepos()
