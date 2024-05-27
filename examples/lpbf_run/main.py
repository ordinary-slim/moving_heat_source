import gmsh
import numpy as np
import MovingHeatSource as mhs
from MovingHeatSource.adaptiveStepper import AdaptiveStepper

inputFile = "input.yaml"
problemInput = mhs.readInput( inputFile )

substrateZLen = 1
partZLen = 1
xLength = 5.0
yLength = 5.0
halfLenX = xLength / 2
halfLenY = yLength / 2
laserRadius = problemInput["radius"]
laserDiameter = 2*laserRadius
layerThickness = problemInput["layerThickness"]

def writeGcode(fileName="Path.gcode"):
    '''
    Write serpentine path
    '''
    gcodeLines = []
    nLayers = int(partZLen / layerThickness)
    Z = 0.0
    Y = -halfLenY
    X = -halfLenX
    E = 0.0
    toggleX = True# if False, toggling in Y
    gcodeLines.append( "G0 F10" )
    for ilayer in range(nLayers):
        # Start with a G0 command
        gcodeLines.append( "G0 X{:.2f} Y{:.2f} Z{:.2f}".format(X, Y, Z) )
        nLongHatches = -1
        length = None
        coord = None
        delta = None
        if ilayer%2==0:
            length = xLength
            coord = Y
            toggleX = True
        else:
            length = yLength
            coord = X
            toggleX = False
        if coord < 0:
            delta = +laserRadius
        else:
            delta = -laserRadius
        nLongHatches = length / laserRadius + 1
        nLongHatches = int(nLongHatches)
        for ihatch in range( nLongHatches ):
            # LONG HATCH
            # Toggle
            if toggleX:
                X *= -1
            else:
                Y *= -1
            E += 0.1
            gcodeLines.append( "G1 X{:.2f} Y{:.2f} E{:.2f}".format(X, Y, E) )
            if ihatch+1 < nLongHatches:
                # SHORT HATCH
                if toggleX:
                    Y += delta
                else:
                    X += delta
                E += 0.02
                gcodeLines.append( "G0 X{:.2f} Y{:.2f} E{:.2f}".format(X, Y, E) )
        Z += layerThickness

    with open(fileName, "w") as gcodeFile:
        gcodeFile.writelines( [line+"\n" for line in gcodeLines] )

    return fileName



def getGmshModel(xLength, yLength, partZLen, substrateZLen):
    gmsh.model.add("prism")
    zLength = substrateZLen + partZLen
    botFaceCenter = np.array([0.0, 0.0, zLength/2])
    gmsh.model.geo.addPoint(botFaceCenter[0]-xLength/2,botFaceCenter[0]-yLength/2,-substrateZLen)
    gmsh.model.geo.addPoint(botFaceCenter[0]+xLength/2,botFaceCenter[0]-yLength/2,-substrateZLen)
    gmsh.model.geo.addPoint(botFaceCenter[0]+xLength/2,botFaceCenter[0]+yLength/2,-substrateZLen)
    gmsh.model.geo.addPoint(botFaceCenter[0]-xLength/2,botFaceCenter[0]+yLength/2,-substrateZLen)
    gmsh.model.geo.addLine(4,3)
    gmsh.model.geo.addLine(3,2)
    gmsh.model.geo.addLine(2,1)
    gmsh.model.geo.addLine(1,4)
    gmsh.model.geo.addCurveLoop([2,3,4,1])
    gmsh.model.geo.addPlaneSurface([1])
    extrudedEntities = gmsh.model.geo.extrude( [(2, 1)], 0.0, 0.0, zLength )

    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(3, [1], name="Domain")

    # TRANSFINITE
    # CURVES
    ## BOT
    nelsInPlaneLines = 50
    nelsVerticalLines = 20
    gmsh.model.geo.mesh.setTransfiniteCurve(1, nelsInPlaneLines+1)
    gmsh.model.geo.mesh.setTransfiniteCurve(2, nelsInPlaneLines+1)
    gmsh.model.geo.mesh.setTransfiniteCurve(3, nelsInPlaneLines+1)
    gmsh.model.geo.mesh.setTransfiniteCurve(4, nelsInPlaneLines+1)
    ## TOP
    gmsh.model.geo.mesh.setTransfiniteCurve(6, nelsInPlaneLines+1)
    gmsh.model.geo.mesh.setTransfiniteCurve(7, nelsInPlaneLines+1)
    gmsh.model.geo.mesh.setTransfiniteCurve(8, nelsInPlaneLines+1)
    gmsh.model.geo.mesh.setTransfiniteCurve(9, nelsInPlaneLines+1)
    ## VERTICAL
    gmsh.model.geo.mesh.setTransfiniteCurve(11, nelsVerticalLines+1)
    gmsh.model.geo.mesh.setTransfiniteCurve(12, nelsVerticalLines+1)
    gmsh.model.geo.mesh.setTransfiniteCurve(16, nelsVerticalLines+1)
    gmsh.model.geo.mesh.setTransfiniteCurve(20, nelsVerticalLines+1)
    # SURFACES
    gmsh.model.geo.mesh.setTransfiniteSurface( 1 )
    for dim, tag in extrudedEntities:
        if dim==2:
            gmsh.model.geo.mesh.setTransfiniteSurface( tag )
        elif dim==3:
            gmsh.model.geo.mesh.setTransfiniteVolume( tag )

    # Hexahedral els
    gmsh.option.setNumber("Mesh.RecombineAll", 1)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    return gmsh.model

def deactivateBelowSurface(p, surfaceZ = 0):
    nels = p.domain.mesh.nels
    activeEls = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElement( ielem )
        if (e.getCentroid()[2] < surfaceZ):
            activeEls.append( ielem )
    substrateEls = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, activeEls )
    p.domain.setActivation( substrateEls )


if __name__=="__main__":
    # MESH
    gmsh.initialize()
    gmshModel = getGmshModel(xLength, yLength, partZLen, substrateZLen)
    meshInput = {}
    meshInput["points"], meshInput["cells"], meshInput["cell_type"] \
        = mhs.gmshModelToMesh( gmshModel )
    gmsh.finalize()
    mesh = mhs.Mesh( meshInput )

    # PROBLEM INPUT

    # read input
    pFixed = mhs.Problem(mesh, problemInput, caseName="fixed")
    deactivateBelowSurface( pFixed, surfaceZ=0 )

    # Set laser path
    gcodeFile = writeGcode(problemInput["path"])

    myDriver = AdaptiveStepper( pFixed,
                                      adimMaxSubdomainSize=10,
                                      threshold=0.1,
                                      meshDenMoving=20,
                                      )

    while not(pFixed.mhs.path.isOver(pFixed.time)):
        myDriver.iterate()
