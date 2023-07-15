import MovingHeatSource as mhs
import numpy as np
import meshzoo
import pdb

def mesh(box, meshDen=4):
    cell_type="quad4"
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

def setAdimR( adimR, input ):
    r = input["radius"]
    HeatSourceSpeedX = max( abs(input["HeatSourceSpeedX"]), abs(input["advectionSpeedX"]))
    HeatSourceSpeedY = max( abs(input["HeatSourceSpeedY"]), abs(input["advectionSpeedY"]))
    HeatSourceSpeedZ = max( abs(input["HeatSourceSpeedZ"]), abs(input["advectionSpeedZ"]))
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
    inputFile = "input.yaml"
    boxPhys = [-16, 16, -5, 5]
    boxBg = [-32, 32, -5, 5]
    adimR = 1

    # read input
    problemInput = mhs.readInput( inputFile )

    FineFRFInput = dict( problemInput )
    FRFInput = dict( problemInput )
    NoTransportMRFInput = dict( problemInput )
    TransportedMRFInput = dict( problemInput )

    # Mesh
    meshInputPhys, meshInputBg = {}, {}
    meshInputPhys["points"], meshInputPhys["cells"], meshInputPhys["cell_type"] = mesh(boxPhys)
    meshInputBg["points"], meshInputBg["cells"], meshInputBg["cell_type"] = mesh(boxBg)

    meshFineFRF       = mhs.Mesh(meshInputPhys)
    meshFRF           = mhs.Mesh(meshInputPhys)
    meshNoTransportMRF= mhs.Mesh(meshInputBg)
    meshTransportedMRF= mhs.Mesh(meshInputBg)
    meshMRFTransporter= mhs.Mesh(meshInputPhys)

    # mhs.Problem params
    # set dt
    dt = setAdimR( adimR, FRFInput )
    for input in [FRFInput, TransportedMRFInput, NoTransportMRFInput]:
        input["dt"] = dt
    ##determine fine problem tstep size
    approxFine_dt = pow(dt, 2)
    approxFine_dt = min( approxFine_dt, dt / 32.0 )
    fineStepsPerStep = int( np.ceil( dt / approxFine_dt ) )
    fine_dt = dt / float( fineStepsPerStep )
    FineFRFInput["dt"] = fine_dt

    #set MRF business NO TRANSPORT
    NoTransportMRFInput["isAdvection"] = 1
    NoTransportMRFInput["advectionSpeedX"] = -FRFInput["HeatSourceSpeedX"]
    NoTransportMRFInput["speedFRF_X"]      = FRFInput["HeatSourceSpeedX"]
    NoTransportMRFInput["HeatSourceSpeedX"] = 0.0
    #set MRF business TRANSPORT
    TransportedMRFInput["isAdvection"] = 1
    TransportedMRFInput["advectionSpeedX"] = -FRFInput["HeatSourceSpeedX"]
    TransportedMRFInput["speedFRF_X"]      = FRFInput["HeatSourceSpeedX"]
    TransportedMRFInput["HeatSourceSpeedX"] = 0.0

    pFineFRF         = mhs.Problem(meshFineFRF, FineFRFInput, caseName="fineFRF")
    pFRF             = mhs.Problem(meshFRF, FRFInput, caseName="FRF")
    pNoTransportMRF  = mhs.Problem(meshNoTransportMRF, NoTransportMRFInput, caseName="NoTransportMRF")
    pTransportedMRF  = mhs.Problem(meshTransportedMRF, TransportedMRFInput, caseName="TransportedMRF")
    pMRFTransporter  = mhs.Problem(meshMRFTransporter, FRFInput, caseName="MRFTransporter")

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
        activeElements = isInsideBox( pNoTransportMRF.domain.mesh, boxPhys )
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

        activeElements = isInsideBox( pTransportedMRF.domain.mesh, boxPhys )
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
