import MovingHeatSource as mhs
import numpy as np
import meshzoo

speed = 10.0

def mesh(box, meshDen=1):
    cell_type="quad4"
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
        np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
        cell_type=cell_type,
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( np.uint32 )
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

def run():
    inputFile = "input.py"
    boxPhys = [-16, 16, -5, 5]
    boxBg = [-32, 32, -5, 5]

    # read input
    problemInput = mhs.readInput( inputFile )

    FRFInput = dict( problemInput )
    TransportedMRFInput = dict( problemInput )

    # Mesh
    meshInputPhys, meshInputBg = {}, {}
    meshInputPhys["points"], meshInputPhys["cells"], meshInputPhys["cell_type"] = mesh(boxPhys)
    meshInputBg["points"], meshInputBg["cells"], meshInputBg["cell_type"] = mesh(boxBg)

    meshTransportedMRF= mhs.Mesh(meshInputBg)
    meshMRFTransporter= mhs.Mesh(meshInputPhys)

    # mhs.Problem params
    # set dt
    dt = 0.5
    TransportedMRFInput["dt"] = dt

    #set MRF business NO TRANSPORT
    TransportedMRFInput["isAdvection"] = 1
    TransportedMRFInput["advectionSpeedX"] = -speed
    TransportedMRFInput["speedFRF_X"]      = speed
    TransportedMRFInput["HeatSourceSpeedX"] = 0.0

    pTransportedMRF  = mhs.Problem(meshTransportedMRF, TransportedMRFInput, caseName="TransportedMRF")
    pMRFTransporter  = mhs.Problem(meshMRFTransporter, FRFInput, caseName="MRFTransporter")

    pMRFTransporter.unknown.interpolate(  pTransportedMRF.unknown )
    for tF, sF in zip(pMRFTransporter.previousValues, pTransportedMRF.previousValues):
        tF.interpolate( sF )

    maxIter = 2
    # FORWARD
    for iteration in range(maxIter):
        #iter transportedMRF
        pTransportedMRF.preiterate(False)#motion
        activeElements = isInsideBox( pTransportedMRF.domain.mesh, boxPhys )
        pTransportedMRF.domain.setActivation( activeElements )

        pTransportedMRF.unknown.interpolate( pMRFTransporter.unknown )
        for tF, sF in zip(pTransportedMRF.previousValues, pMRFTransporter.previousValues):
            tF.interpolate( sF )


        pTransportedMRF.preAssemble(False)#update forced dofs, reallocate
        pTransportedMRF.iterate()

        pMRFTransporter.fakeIter()
        pMRFTransporter.unknown.interpolate( pTransportedMRF.unknown )
        for tF, sF in zip(pMRFTransporter.previousValues, pTransportedMRF.previousValues):
            tF.interpolate( sF )
        pTransportedMRF.writepos()

def test():
    run()
    referenceFile = "post_TransportedMRF_reference/TransportedMRF_2.vtu"
    newFile = "post_TransportedMRF/TransportedMRF_2.vtu"
    assert mhs.meshio_comparison(referenceFile, newFile)

if __name__=="__main__":
    test()
