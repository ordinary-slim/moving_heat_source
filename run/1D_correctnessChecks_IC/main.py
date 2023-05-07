import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Release/')
import MovingHeatSource as mhs
from wrapper import Problem
import numpy as np
import meshzoo
import pdb
import meshio

speed = 10
bumper = 20
Tfinal = 2.0

def mesh(leftEnd, rightEnd, elDen):
    cell_type="line2"
    nels = int((rightEnd-leftEnd)*elDen)
    points = np.linspace( leftEnd, rightEnd, nels+1, dtype=np.float64)
    points = points.reshape( (nels+1, 1) )
    auxCells = []
    for iel in range(nels):
        auxCells.append( [iel, iel+1] )
    cells = np.array( auxCells, dtype=int )
    return points, cells, cell_type

def isInside( mesh, box ):
    activeElements = []
    for ielem in range( mesh.nels ):
        el = mesh.getElement( ielem )
        pos = mesh.posFRF[el.con]
        xmin = min(pos[:, 0])
        xmax = max(pos[:, 0])

        isInside = 1*(xmin>=box[0] and xmax <= box[1])
        activeElements.append(isInside)
    return activeElements

if __name__=="__main__":
    inputFile = "input.txt"

    pFRF = Problem("FRF")
    pMRF_Advec = Problem("MRF_Advec")
    pMRF_Trans = Problem("MRF_Trans")
    pMRF_TransHelper = Problem("TransHelper")
    pMRF_AdvecTrans = Problem("MRF_AdvecTrans")
    pMRF_AdvecTransHelper = Problem("AdvecTransHelper")
    pFRF_Step  = Problem("FRFStep")

    elDen = 16
    dt = 0.5
    leftEnd = -25.0
    rightEnd = +25.0
    L = rightEnd - leftEnd
    boxDomain = [leftEnd, rightEnd]
    boxBG = [leftEnd-bumper, rightEnd+bumper]
    # FRF meshes
    points, cells, cell_type = mesh(boxDomain[0], boxDomain[1], elDen)
    for p in [pFRF, pMRF_TransHelper, pMRF_AdvecTransHelper, pFRF_Step]:
        p.input["points"] = points
        p.input["cells"] = cells
        p.input["cell_type"]=cell_type
        #p.input["numberOfGaussPointsCells"]=3

    # background meshes
    points, cells, cell_type = mesh(boxBG[0], boxBG[1], elDen)
    for p in [pMRF_Advec, pMRF_Trans, pMRF_AdvecTrans]:
        p.input["points"] = points
        p.input["cells"] = cells
        p.input["cell_type"]=cell_type
        #p.input["numberOfGaussPointsCells"]=3

    #read input
    for p in [pFRF, pMRF_Advec, pMRF_Trans, pMRF_TransHelper, pMRF_AdvecTrans, pMRF_AdvecTransHelper, pFRF_Step]:
        p.parseInput( inputFile )
        p.input["dt"] = dt

    #set MRF business TRANSPORT
    for p in [pMRF_Trans, pMRF_AdvecTrans]:
        p.input["speedFRF_X"]  = speed
        p.input["HeatSourceSpeedX"] = -speed
        p.input["isStabilized"] = 0
    #set MRF business ADVEC
    for p in [pMRF_Advec, pMRF_AdvecTrans, pMRF_AdvecTrans]:
        p.input["speedFRF_X"]  = speed
        p.input["HeatSourceSpeedX"] = -speed
        p.input["advectionSpeedX"] = -speed
        p.input["isStabilized"] = 1

    for p in [pFRF, pMRF_Advec, pMRF_Trans, pMRF_TransHelper, pMRF_AdvecTrans, pMRF_AdvecTransHelper, pFRF_Step]:
        p.initialize()

    # Initial condition
    # Different IC
    #f = lambda pos : abs(pos[0]+pos[1])
    #f = lambda pos : 25 + (1 / (1 +pos[0]**2) )
    #f = lambda pos : max(25, 26-abs(pos[0])/10)
    # Sine with noise
    A, B = 2, 0.05
    freq = 1
    np.random.seed(1)
    f = lambda pos : A*np.sin(freq*2*np.pi/L*pos[0]) + B*np.random.normal(scale=8, size=pos[0].size)
    for p in [pFRF, pMRF_Advec, pMRF_Trans, pMRF_TransHelper, pMRF_AdvecTrans, pMRF_AdvecTransHelper, pFRF_Step]:
        p.forceState( f )
        p.writepos()
    #quickfix

    # FORWARD
    postMRF = True
    shift = None
    while (pFRF.time < Tfinal):
        #Advected MRF
        pMRF_Advec.updateFRFpos()
        pMRF_Advec.iterate()
        pMRF_Advec.writepos()

        if postMRF:
            shift=(-pMRF_Advec.mesh.shiftFRF)

        #FRF
        pFRF.iterate()#assembly + solve
        pFRF.writepos(shift=shift)
        
        # Copy FRF into MRF using the copy constructor
        pMRF_Step = Problem(  "MRFStep", problem=pFRF_Step )
        # Change reference frame
        pMRF_Step.frf2mrf(speed=np.array([speed, 0.0, 0.0], dtype=np.float64))
        pMRF_Step.setStabilization(True)
        pMRF_Step.updateFRFpos()
        pMRF_Step.iterate()
        pMRF_Step.writepos(shift=None)

        # Force FRF and solve
        pFRF_Step.unknown.interpolate2dirichlet( pMRF_Step.unknown )
        pFRF_Step.iterate()
        pFRF_Step.unknown.releaseDirichlet()
        pFRF_Step.writepos(shift=None)
