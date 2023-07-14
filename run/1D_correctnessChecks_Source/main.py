import MovingHeatSource as mhs
import numpy as np
import meshzoo
import pdb
import meshio

speed = 10
bumper = 20
Tfinal = bumper/speed/2

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
    inputFile = "input.py"

    pFRF = mhs.Problem("FRF")
    pMRF_Advec = mhs.Problem("MRF_Advec")
    pMRF_Trans = mhs.Problem("MRF_Trans")
    pMRF_TransHelper = mhs.Problem("TransHelper")
    pMRF_AdvecTrans = mhs.Problem("MRF_AdvecTrans")
    pMRF_AdvecTransHelper = mhs.Problem("AdvecTransHelper")

    elDen = 10
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
        #p.input["numberOfGaussPointsCells"]=3

    # MRF meshes
    points, cells, cell_type = mesh(boxBG[0], boxBG[1], elDen)
    for p in [pMRF_Advec, pMRF_Trans, pMRF_AdvecTrans]:
        p.input["points"] = points
        p.input["cells"] = cells
        p.input["cell_type"]=cell_type
        #p.input["numberOfGaussPointsCells"]=3

    #read input
    for p in [pFRF, pMRF_Advec, pMRF_Trans, pMRF_TransHelper, pMRF_AdvecTrans, pMRF_AdvecTransHelper]:
        p.parseInput( inputFile )

    print( "h = {} , dx = {}, dt = {}".format( points[1]-points[0], speed * pFRF.input["dt"], pFRF.input["dt"]) )

    #set MRF business ADVEC
    for p in [pMRF_Advec, pMRF_AdvecTrans, pMRF_AdvecTrans]:
        p.input["speedFRF_X"]  = speed
        p.input["HeatSourceSpeedX"] = -speed
        p.input["advectionSpeedX"] = -speed
        p.input["isStabilized"] = 0
    #set MRF business TRANSPORT
    for p in [pMRF_Trans, pMRF_AdvecTrans]:
        p.input["speedFRF_X"]  = speed
        p.input["HeatSourceSpeedX"] = -speed

    for p in [pFRF, pMRF_Advec, pMRF_Trans, pMRF_TransHelper, pMRF_AdvecTrans, pMRF_AdvecTransHelper]:
        p.initialize()

    pMRF_TransHelper.unknown.interpolate(  pMRF_Trans.unknown )
    pMRF_AdvecTransHelper.unknown.interpolate(  pMRF_AdvecTrans.unknown )

    # FORWARD
    while (pFRF.time < Tfinal):
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
        pMRF_Trans.unknown.interpolate( pMRF_TransHelper.unknown )
        activeElements = isInside( pMRF_Trans.mesh, [leftEnd, rightEnd])
        pMRF_Trans.activate( activeElements )
        pMRF_Trans.iterate()
        pMRF_TransHelper.fakeIter()
        pMRF_TransHelper.unknown.interpolate( pMRF_Trans.unknown )
        pMRF_Trans.writepos()
        pMRF_TransHelper.writepos()

        #Advection + transport
        pMRF_AdvecTrans.updateFRFpos()
        pMRF_AdvecTrans.unknown.interpolate( pMRF_AdvecTransHelper.unknown )
        activeElements = isInside( pMRF_AdvecTrans.mesh, [leftEnd, rightEnd])
        pMRF_AdvecTrans.activate( activeElements )
        pMRF_AdvecTrans.iterate()
        pMRF_AdvecTransHelper.fakeIter()
        pMRF_AdvecTransHelper.unknown.interpolate( pMRF_AdvecTrans.unknown )
        pMRF_AdvecTrans.writepos()
        pMRF_AdvecTransHelper.writepos()

