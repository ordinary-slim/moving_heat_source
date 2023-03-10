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
Tfinal = bumper/speed/2
maxIter = 100

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

    elDen = 20
    leftEnd = -50.0
    rightEnd = +50.0
    boxDomain = [leftEnd, rightEnd]
    boxBG = [leftEnd-bumper, rightEnd+bumper]
    # FRF meshes
    points, cells, cell_type = mesh(boxDomain[0], boxDomain[1], elDen)
    for p in [pFRF, pMRF_TransHelper, pMRF_AdvecTransHelper]:
        p.input["points"] = points
        p.input["cells"] = cells
        p.input["cell_type"]=cell_type
        #p.input["numberOfGaussPoints"]=3

    # MRF meshes
    points, cells, cell_type = mesh(boxBG[0], boxBG[1], elDen)
    for p in [pMRF_Advec, pMRF_Trans, pMRF_AdvecTrans]:
        p.input["points"] = points
        p.input["cells"] = cells
        p.input["cell_type"]=cell_type
        #p.input["numberOfGaussPoints"]=3

    #read input
    for p in [pFRF, pMRF_Advec, pMRF_Trans, pMRF_TransHelper, pMRF_AdvecTrans, pMRF_AdvecTransHelper]:
        p.parseInput( inputFile )

    #set MRF business TRANSPORT
    for p in [pMRF_Trans, pMRF_AdvecTrans]:
        p.input["speedFRF_X"]  = speed
        p.input["speedX"] = -speed
        p.input["isStabilized"] = 0
    #set MRF business ADVEC
    for p in [pMRF_Advec, pMRF_AdvecTrans, pMRF_AdvecTrans]:
        p.input["speedFRF_X"]  = speed
        p.input["speedX"] = -speed
        p.input["advectionSpeedX"] = -speed
        p.input["isStabilized"] = 1

    for p in [pFRF, pMRF_Advec, pMRF_Trans, pMRF_TransHelper, pMRF_AdvecTrans, pMRF_AdvecTransHelper]:
        p.initialize()

    # Initial condition
    # Different IC
    #f = lambda pos : abs(pos[0]+pos[1])
    f = lambda pos : 25 + (1 / (1 +pos[0]**2) )
    f = lambda pos : max(25, 26-abs(pos[0])/10)
    for p in [pFRF, pMRF_Advec, pMRF_Trans, pMRF_TransHelper, pMRF_AdvecTrans, pMRF_AdvecTransHelper]:
        p.forceState( f )
    #quickfix

    # FORWARD
    iteration = 0
    while ((pFRF.time < Tfinal)and(iteration < maxIter)):
        iteration += 1
        #FRF
        pFRF.iterate()#assembly + solve
        pFRF.writepos()
        
        #Advected MRF
        pMRF_Advec.updateFRFpos()
        #pdb.set_trace()
        pMRF_Advec.iterate()
        pMRF_Advec.writepos()

        #Transported MRF
        pMRF_Trans.updateFRFpos()
        pMRF_Trans.unknown.getFromExternal( pMRF_TransHelper.unknown )
        activeElements = isInside( pMRF_Trans.mesh, [leftEnd, rightEnd])
        pMRF_Trans.activate( activeElements )
        pMRF_Trans.iterate()
        pMRF_TransHelper.fakeIter()
        pMRF_TransHelper.unknown.getFromExternal( pMRF_Trans.unknown )
        pMRF_Trans.writepos()
        pMRF_TransHelper.writepos()

        #Advection + transport
        pMRF_AdvecTrans.updateFRFpos()
        pMRF_AdvecTrans.unknown.getFromExternal( pMRF_AdvecTransHelper.unknown )
        activeElements = isInside( pMRF_AdvecTrans.mesh, [leftEnd, rightEnd])
        pMRF_AdvecTrans.activate( activeElements )
        pMRF_AdvecTrans.iterate()
        pMRF_AdvecTransHelper.fakeIter()
        pMRF_AdvecTransHelper.unknown.getFromExternal( pMRF_AdvecTrans.unknown )
        pMRF_AdvecTrans.writepos()
        pMRF_AdvecTransHelper.writepos()
