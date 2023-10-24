import meshzoo
import numpy as np
import MovingHeatSource as mhs
import gmsh
import pdb

problemInput = mhs.readInput( "input.yaml" )
partLens = np.array(problemInput["part"])
substrateLens = np.array(problemInput["substrate"])
radiusHs = problemInput["radius"]
elsPerRadius = problemInput["elsPerRadius"]
fineElSize = radiusHs / elsPerRadius
nBounLayers = 4
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

def getMeshPhysical(popup=False):
    gmsh.initialize()
    halfLensPart = np.array( partLens ) / 2
    halfLensSubstrate = np.array( substrateLens ) / 2
    # Bot surface part
    gmsh.model.geo.addPoint( -halfLensPart[0], 0.0, 0.0, tag=1 )
    gmsh.model.geo.addPoint( -halfLensPart[0], +halfLensPart[1], 0.0, tag=2 )
    gmsh.model.geo.addPoint( +halfLensPart[0], +halfLensPart[1], 0.0, tag=3 )
    gmsh.model.geo.addPoint( +halfLensPart[0], 0.0, 0.0, tag=4 )
    linesBottomSurfacePart = []
    linesBottomSurfacePart.append( gmsh.model.geo.addLine( 1, 2, tag=1 ) )
    linesBottomSurfacePart.append( gmsh.model.geo.addLine( 2, 3, tag=2 ) )
    linesBottomSurfacePart.append( gmsh.model.geo.addLine( 3, 4, tag=3 ) )
    linesBottomSurfacePart.append( gmsh.model.geo.addLine( 4, 1, tag=4 ) )
    for idx, line in enumerate(linesBottomSurfacePart):
        numEls = halfLensPart[1] / fineElSize
        if idx%2 != 0:
            numEls = partLens[0] / fineElSize
        gmsh.model.geo.mesh.setTransfiniteCurve(line, round(numEls)+1 )
    curveLoopBotSurfacePart = gmsh.model.geo.addCurveLoop( linesBottomSurfacePart, 1 )
    botSurfacePart = gmsh.model.geo.addPlaneSurface([curveLoopBotSurfacePart], 1)
    gmsh.model.geo.mesh.setTransfiniteSurface(botSurfacePart)
    gmsh.model.geo.mesh.setRecombine(2,botSurfacePart) 

    # Substrate extrusions
    ## Substrate extrusions Z
    ## Uniform extrusion
    lenUniformBotExtrusion = halfLensPart[2]
    numElements = round(lenUniformBotExtrusion / fineElSize )
    uniformBotExtrusion = gmsh.model.geo.extrude([(2, botSurfacePart)], 0, 0, -lenUniformBotExtrusion, numElements=[numElements], recombine=True)
    midBotSurface, _, xSideMidBot1, ySideMidBot1, xSideMidBot2, ySideMidBot2 = uniformBotExtrusion
    ## Coarsening
    lenCoarseBotExtrusion = substrateLens[2] - nBounLayers*fineElSize
    nElements, heights = boundaryLayerProgression( lenCoarseBotExtrusion, nBounLayers, fineElSize, coarseElFactor=coarseElFactor )
    coarseBotExtrusion = gmsh.model.geo.extrude([(2, midBotSurface[1])], 0, 0, -lenCoarseBotExtrusion, numElements=nElements, heights=heights, recombine=True)
    botBotSurface, _, xSideBotBot1, ySideBotBot1, xSideBotBot2, ySideBotBot2 = coarseBotExtrusion

    ## Substrate extrusions Y
    ### Uniform extrusions
    '''
    uniformYExtrusions = []
    lenUniformYExtrusion = nBounLayers*fineElSize
    numElements = np.round(lenUniformYExtrusion / fineElSize).astype(int)
    for idx, surface in enumerate([ySideMidBot1, ySideBotBot1]):
        uniformYExtrusions.append( gmsh.model.geo.extrude([(2, surface[1])], 0.0, lenUniformYExtrusion, 0.0, numElements=[numElements], recombine=True) )
    '''

    ### Coarse extrusions
    coarseYExtrusions = []
    lenCoarseYExtrusion = halfLensSubstrate[1] - halfLensPart[1]
    nElements, heights = boundaryLayerProgression( lenCoarseYExtrusion, nBounLayers, fineElSize, coarseElFactor=coarseElFactor )
    for idx, surface in enumerate([ySideMidBot1, ySideBotBot1]):
        coarseYExtrusions.append( gmsh.model.geo.extrude([(2, surface[1])], 0.0, lenCoarseYExtrusion, 0.0, numElements=nElements, heights=heights, recombine=True) )

    ## Substrate extrusions X
    ### Uniform extrusions
    positiveXExtrusions = []
    negativeXExtrusions = []
    surfacesXPlus = []
    surfacesXMinus = []
    # Collect all X surfaces
    for extrusion in [uniformBotExtrusion, coarseBotExtrusion]:
        surfacesXMinus.append( extrusion[2][1] )
        surfacesXPlus.append( extrusion[4][1] )
    for idx, extrusion in enumerate(coarseYExtrusions):
        surfacesXPlus.append( extrusion[3][1] )
        surfacesXMinus.append( extrusion[5][1] )
    extrusionLen = (substrateLens[0] - partLens[0]) / 2
    nElements, heights = boundaryLayerProgression( extrusionLen, nBounLayers, fineElSize, coarseElFactor=coarseElFactor )
    for surface in surfacesXPlus:
        positiveXExtrusions.append( gmsh.model.geo.extrude([(2, surface)], extrusionLen, 0.0, 0.0, numElements=nElements, heights=heights, recombine=True ) )
    for surface in surfacesXMinus:
        negativeXExtrusions.append( gmsh.model.geo.extrude([(2, surface)], -extrusionLen, 0.0, 0.0, numElements=nElements, heights=heights, recombine=True ) )
    
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
