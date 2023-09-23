import meshzoo
import numpy as np
import MovingHeatSource as mhs
import gmsh

problemInput = mhs.readInput( "input.yaml" )
partLens = problemInput["part"]
substrateLens = problemInput["substrate"]
radiusHs = problemInput["radius"]
layerThickness = problemInput["layerThickness"]
fineElSize = min(radiusHs, layerThickness)
nelsPerLayer = 1
nBounLayers = 2
coarseElFactor = 4

def boundaryLayerProgression( extrusionSize, nBounLayers, fineElSize, coarseElFactor=2 ):
    approxCoarseElSize = coarseElFactor * fineElSize
    numCoarseEls = np.ceil( (extrusionSize - nBounLayers*fineElSize) / approxCoarseElSize )
    coarsestElSize = (extrusionSize - nBounLayers*fineElSize) / numCoarseEls
    numMidCoarseEls = np.ceil( numCoarseEls / 2 )
    midCoarseElSize = coarsestElSize / 2.0
    numCoarsestEls = numCoarseEls - numMidCoarseEls

    numElementsPerLayer = np.array( [nBounLayers] + [numMidCoarseEls] + [numCoarsestEls])
    heights = np.array( [fineElSize, midCoarseElSize, coarsestElSize] )

    # Format heights to cumsum
    heights *= numElementsPerLayer
    heights = np.cumsum( heights )
    heights /= heights[-1]

    return numElementsPerLayer, heights

def meshBox(box, elementSize=0.25):
    cell_type="hexa8"
    nelsX = int((box[1] - box[0])/elementSize)
    nelsY = int((box[3] - box[2])/elementSize)
    nelsZ = int((box[5] - box[4])/elementSize)

    points, cells = meshzoo.cube_hexa(
        np.linspace( box[0], box[1], nelsX+1),
        np.linspace( box[3], box[2], nelsY+1),
        np.linspace( box[5], box[4], nelsZ+1),
    )
    cells = cells.astype( np.uint32 )
    meshDict = {}
    meshDict["points"], meshDict["cells"], meshDict["cell_type"] = points, cells, cell_type
    meshDict["dimension"] = 3
    return mhs.Mesh( meshDict )

def meshAroundHS( adimR, problemInput, elementSize=0.25 ):
    radius = problemInput["radius"]
    initialPosition = problemInput["initialPosition"]
    trailLength = adimR * radius
    capotLength = min( trailLength, 2*radius )
    halfLengthY = min( trailLength, capotLength )
    halfLengthZ = halfLengthY
    box = [initialPosition[0] - trailLength, initialPosition[0] + capotLength,
           initialPosition[1] - halfLengthY, initialPosition[1] + halfLengthY,
           initialPosition[2] - halfLengthZ, initialPosition[2] + halfLengthZ,
           ]

    return meshBox( box, elementSize = elementSize )

def getMeshPhysicaMeshzoo():
    return meshBox( [-25, 25, -5, 5, -5, 1], elementSize=1 )

def getMeshPhysical(popup=False):
    gmsh.initialize()
    halfLensPart = np.array( partLens ) / 2
    halfLensSubstrate = np.array( substrateLens ) / 2
    # Bot surface part
    gmsh.model.geo.addPoint( -halfLensPart[0], -halfLensPart[1], 0.0, tag=1 )
    gmsh.model.geo.addPoint( -halfLensPart[0], +halfLensPart[1], 0.0, tag=2 )
    gmsh.model.geo.addPoint( +halfLensPart[0], +halfLensPart[1], 0.0, tag=3 )
    gmsh.model.geo.addPoint( +halfLensPart[0], -halfLensPart[1], 0.0, tag=4 )
    linesBottomSurfacePart = []
    linesBottomSurfacePart.append( gmsh.model.geo.addLine( 1, 2, tag=1 ) )
    linesBottomSurfacePart.append( gmsh.model.geo.addLine( 2, 3, tag=2 ) )
    linesBottomSurfacePart.append( gmsh.model.geo.addLine( 3, 4, tag=3 ) )
    linesBottomSurfacePart.append( gmsh.model.geo.addLine( 4, 1, tag=4 ) )
    for idx, line in enumerate(linesBottomSurfacePart):
        numEls = partLens[1] / fineElSize
        if idx%2 != 0:
            numEls = partLens[0] / fineElSize
        gmsh.model.geo.mesh.setTransfiniteCurve(line, round(numEls)+1 )
    curveLoopBotSurfacePart = gmsh.model.geo.addCurveLoop( linesBottomSurfacePart, 1 )
    botSurfacePart = gmsh.model.geo.addPlaneSurface([curveLoopBotSurfacePart], 1)
    gmsh.model.geo.mesh.setTransfiniteSurface(botSurfacePart)
    gmsh.model.geo.mesh.setRecombine(2,botSurfacePart) 
    # Part extrusion
    nelsPartZ = partLens[2] / (layerThickness / nelsPerLayer) # 2 el / layer
    topExtrusion = gmsh.model.geo.extrude([(2, botSurfacePart)], 0, 0, +partLens[2], numElements=[nelsPartZ], recombine=True)
    # Substrate extrusions
    ## Substrate extrusions Z
    ## Uniform extrusion
    lenUniformBotExtrusion = nBounLayers*fineElSize
    uniformBotExtrusion = gmsh.model.geo.extrude([(2, botSurfacePart)], 0, 0, -lenUniformBotExtrusion, numElements=[nBounLayers], recombine=True)
    midBotSurface, _, xSideMidBot1, ySideMidBot1, xSideMidBot2, ySideMidBot2 = uniformBotExtrusion
    ## Coarsening
    lenCoarseBotExtrusion = substrateLens[2] - nBounLayers*fineElSize
    nElements, heights = boundaryLayerProgression( lenCoarseBotExtrusion, nBounLayers, fineElSize, coarseElFactor=coarseElFactor )
    coarseBotExtrusion = gmsh.model.geo.extrude([(2, midBotSurface[1])], 0, 0, -lenCoarseBotExtrusion, numElements=nElements, heights=heights, recombine=True)
    botBotSurface, _, xSideBotBot1, ySideBotBot1, xSideBotBot2, ySideBotBot2 = coarseBotExtrusion
    ## Substrate extrusions Y
    ### Uniform extrusions
    uniformYExtrusions = []
    lenUniformYExtrusion = nBounLayers*fineElSize
    for idx, surface in enumerate([ySideMidBot1, ySideMidBot2, ySideBotBot1, ySideBotBot2]):
        uniformYExtrusions.append( gmsh.model.geo.extrude([(2, surface[1])], 0.0, np.power(-1, idx)*lenUniformYExtrusion, 0.0, numElements=[nBounLayers], recombine=True) )

    ### Coarse extrusions
    coarseYExtrusions = []
    lenCoarseYExtrusion = (substrateLens[1] - partLens[1])/2 - nBounLayers*fineElSize
    nElements, heights = boundaryLayerProgression( lenCoarseYExtrusion, nBounLayers, fineElSize, coarseElFactor=coarseElFactor )
    for idx, extrusion in enumerate( uniformYExtrusions ):
        tagSurface = extrusion[0][1]
        coarseYExtrusions.append( gmsh.model.geo.extrude([(2, tagSurface)], 0.0, np.power(-1, idx)*lenCoarseYExtrusion, 0.0, numElements=nElements, heights=heights, recombine=True) )

    ## Substrate extrusions X
    ### Uniform extrusions
    positiveUniformExtrusionsX = []
    negativeUniformExtrusionsX = []
    surfacesXPlus = []
    surfacesXMinus = []
    # Collect all X surfaces
    for extrusion in [uniformBotExtrusion, coarseBotExtrusion]:
        surfacesXMinus.append( extrusion[2][1] )
        surfacesXPlus.append( extrusion[4][1] )
    for idx, extrusion in enumerate(uniformYExtrusions + coarseYExtrusions):
        if (idx%2 == 0):
            surfacesXPlus.append( extrusion[3][1] )
            surfacesXMinus.append( extrusion[5][1] )
        else:
            surfacesXMinus.append( extrusion[3][1] )
            surfacesXPlus.append( extrusion[5][1] )
    extrusionLen = nBounLayers*fineElSize
    numElements = extrusionLen / fineElSize
    for surface in surfacesXPlus:
        positiveUniformExtrusionsX.append( gmsh.model.geo.extrude([(2, surface)], extrusionLen, 0.0, 0.0, numElements=[numElements], recombine=True ) )
    for surface in surfacesXMinus:
        negativeUniformExtrusionsX.append( gmsh.model.geo.extrude([(2, surface)], -extrusionLen, 0.0, 0.0, numElements=[numElements], recombine=True ) )
    ### Coarse extrusions
    extrusionLen = (substrateLens[0] - partLens[0])/2 - nBounLayers*fineElSize
    nElements, heights = boundaryLayerProgression( extrusionLen, nBounLayers, fineElSize, coarseElFactor=coarseElFactor )
    for extrusion in positiveUniformExtrusionsX:
        gmsh.model.geo.extrude([extrusion[0]], +extrusionLen, 0.0, 0.0, numElements=nElements, heights = heights,  recombine = True )
    for extrusion in negativeUniformExtrusionsX:
        gmsh.model.geo.extrude([extrusion[0]], -extrusionLen, 0.0, 0.0, numElements=nElements, heights = heights, recombine = True )


    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    volumeTags = []
    for _, tag in gmsh.model.getEntities( 3 ):
        volumeTags.append( tag )
    gmsh.model.addPhysicalGroup(3, volumeTags, tag=1, name="Domain")

    if popup:
        gmsh.fltk.run()
    else:
        mesh = mhs.gmshModelToMesh( gmsh.model )
        gmsh.finalize()
        return mesh

if __name__=="__main__":
    getMeshPhysical(popup = True)
