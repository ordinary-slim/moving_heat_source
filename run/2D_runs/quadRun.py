import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
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
    return points, cells

def isInsideBox( mesh, box ):
    activeElements = []
    for ielem in range( mesh.nels ):
        el = mesh.getElement( ielem )
        pos = mesh.pos_noAdv[el.con]
        xmin = min(pos[:, 0])
        xmax = max(pos[:, 0])
        ymin = min(pos[:, 1])
        ymax = max(pos[:, 1])

        isInside = 1*(xmin>=box[0] and xmax <= box[1] and ymin >= box[2] and ymax <= box[3])
        activeElements.append(isInside)
    return activeElements

def computeAvElSize( p ):
    #double counting!
    accumulation = 0
    N = 0
    points = p.input["points"]
    con    = p.input["cells"]
    for ielem in range(con.shape[0]):
        #compute size
        for inode in range( con.shape[1] ):
            accumulation += np.linalg.norm(
                points[con[ielem][(inode+1)%con.shape[1]]] - points[con[ielem][inode]])
    N = con.size
    return accumulation / float( N )

def setAdimR( adimR, p ):
    r = p.input["radius"]
    speedX = max( abs(p.input["speedX"]), abs(p.input["advectionSpeedX"]))
    speedY = max( abs(p.input["speedY"]), abs(p.input["advectionSpeedY"]))
    speedZ = max( abs(p.input["speedZ"]), abs(p.input["advectionSpeedZ"]))
    speed  = np.linalg.norm( np.array( [speedX, speedY, speedZ] ) )
    return (adimR * r / speed)

if __name__=="__main__":
    inputFile = "input.txt"
    boxRef = [-16, 16, -5, 5]
    boxInac = [-36, 36, -5, 5]
    adimR = 1

    problemFine    = Problem("fineFRF")
    problemFRF     = Problem("FRF")
    problemMRF_act = Problem("MRF_act")

    points, cells = mesh(boxRef)
    for p in [problemFine, problemFRF]:
        p.input["points"] = points
        p.input["cells"] = cells
        p.input["cell_type"]="quad4"

    points, cells = mesh(boxInac)
    problemMRF_act.input["points"] = points
    problemMRF_act.input["cells"] = cells
    problemMRF_act.input["cell_type"]="quad4"

    #read input
    for p in [problemFine, problemFRF, problemMRF_act]:
        p.parseInput( inputFile )

    # set dt
    dt = setAdimR( adimR, problemFRF )
    problemFRF.input["dt"] = dt
    problemMRF_act.input["dt"] = dt
    ##determine fine problem tstep size
    approxFine_dt = pow(dt, 2)
    approxFine_dt = min( approxFine_dt, dt / 16.0 )
    fineStepsPerStep = int( np.ceil( dt / approxFine_dt ) )
    fine_dt = dt / float( fineStepsPerStep )
    problemFine.input["dt"] = fine_dt

    #set MRF business
    problemMRF_act.input["isAdvection"] = 1
    problemMRF_act.input["advectionSpeedX"] = -problemMRF_act.input["speedX"]
    problemMRF_act.input["speedX"] = 0.0

    for p in [problemFine, problemFRF, problemMRF_act]:
        p.initialize()

    maxIter = problemFRF.input["maxIter"]
    # FORWARD
    for iteration in range(maxIter):
        #fine problem
        #for istep in range(fineStepsPerStep):
            #problemFine.iterate()
        #problemFine.writepos()

        #for p in [problemFRF, problemMRF_act]:
        for p in [problemFRF, problemMRF_act]:
            p.updateFRFpos()#get tn+1 positions (not tn)
            activeElements = isInsideBox( p.mesh, boxRef )#active tn+1 positions
            p.activate( activeElements )#activation
            p.iterate()#assembly + solve
            p.writepos()

    '''
    problemMRF_act.setAdvectionSpeed( -problemMRF_act.advectionSpeed )
    problemFine.mhs.setSpeed( -problemFine.mhs.speed )
    problemFRF.mhs.setSpeed( -problemFRF.mhs.speed )

    # BACKWARDS
    for iteration in range(maxIter):
        #fine problem
        for istep in range(fineStepsPerStep):
            problemFine.iterate()
        problemFine.writepos()

        #for p in [problemMRF_act]:
        for p in [problemFRF, problemMRF_act]:
            p.updateFRFpos()#get tn+1 positions (not tn)
            activeElements = isInsideBox( p.mesh, boxRef )#active tn+1 positions
            p.activate( activeElements )#activation
            p.iterate()#assembly + solve
            p.writepos()

    '''
