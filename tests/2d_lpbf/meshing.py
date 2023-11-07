import meshzoo
import numpy as np
import MovingHeatSource as mhs
import gmsh
import pdb

problemInput = mhs.readInput( "input.yaml" )
partLens = problemInput["part"]
substrateLens = problemInput["substrate"]
radiusHs = problemInput["radius"]
layerThickness = problemInput["layerThickness"]
fineElSize = min(radiusHs, layerThickness)
nBounLayers = 2

def boundaryLayerProgression( extrusionSize, nBounLayers, fineElSize, coarseElFactor=2 ):
    approxCoarseElSize = coarseElFactor * fineElSize
    numCoarseEls = np.ceil( (extrusionSize - nBounLayers*fineElSize) / approxCoarseElSize )
    coarsestElSize = (extrusionSize - nBounLayers*fineElSize) / numCoarseEls
    numMidCoarseEls = np.ceil( numCoarseEls / 2 )
    midCoarseElSize = coarsestElSize / 2.0
    numCoarsestEls = numCoarseEls - numMidCoarseEls

    numElsPerLayer = [nBounLayers] + [numMidCoarseEls]
    heights = [fineElSize, midCoarseElSize]
    if numCoarsestEls>0:
        numElsPerLayer += [numCoarsestEls]
        heights.append( coarsestElSize ) 
    numElementsPerLayer = np.array(numElsPerLayer)
    heights = np.array(heights)

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
    gmsh.model.geo.addPoint( +halfLensPart[0], 0.0, 0.0, tag=2 )

    lineBottomSurfacePart = gmsh.model.geo.addLine( 1, 2, tag=1 )

    numEls = partLens[0] / fineElSize
    gmsh.model.geo.mesh.setTransfiniteCurve(lineBottomSurfacePart, round(numEls)+1 )

    # Substrate extrusions
    ## Substrate extrusions Z
    ## Uniform extrusion
    lenUniformBotExtrusion = substrateLens[1]
    nelsSubstrateZ = np.round(lenUniformBotExtrusion / fineElSize).astype(int)
    uniformBotExtrusion = gmsh.model.geo.extrude([(1, lineBottomSurfacePart)], 0, 0, -lenUniformBotExtrusion, numElements=[nelsSubstrateZ], recombine=True)
    midBotLine, _, zSideMidBot1, zSideMidBot2 = uniformBotExtrusion

    ## Substrate extrusions X
    ### Uniform extrusions
    uniformXExtrusions = []
    lenUniformYExtrusion = nBounLayers*fineElSize
    for idx, line in enumerate([zSideMidBot1, zSideMidBot2]):
        uniformXExtrusions.append( gmsh.model.geo.extrude([(1, line[1])], np.power(-1, idx)*lenUniformYExtrusion, 0.0, 0.0, numElements=[nBounLayers], recombine=True) )

    ### Coarse extrusions
    coarseXExtrusions = []
    lenCoarseXExtrusion = (substrateLens[0] - partLens[0])/2 - nBounLayers*fineElSize
    nElements, heights = boundaryLayerProgression( lenCoarseXExtrusion, nBounLayers, fineElSize, coarseElFactor=2 )
    for idx, extrusion in enumerate( uniformXExtrusions ):
        tagLine = extrusion[0][1]
        coarseXExtrusions.append( gmsh.model.geo.extrude([(1, tagLine)], np.power(-1, idx)*lenCoarseXExtrusion, 0.0, 0.0, numElements=nElements, heights=heights, recombine=True) )

    # Part extrusion
    topLines = []
    topLines.append( (1, lineBottomSurfacePart) )
    topLines.append( (1, abs(uniformXExtrusions[0][-1][1]) ) )
    topLines.append( (1, abs(uniformXExtrusions[1][-2][1]) ) )
    topLines.append( (1, abs( coarseXExtrusions[0][-1][1]) ) )
    topLines.append( (1, abs( coarseXExtrusions[1][-2][1]) ) )

    nelsPartZ = np.round(partLens[1] / fineElSize).astype(int)
    topExtrusions = []
    for line in topLines:
        topExtrusions.append( gmsh.model.geo.extrude([line], 0, 0, +partLens[1], numElements=[nelsPartZ], recombine=True) )


    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    volumeTags = []
    for _, tag in gmsh.model.getEntities( 2 ):
        volumeTags.append( tag )
    gmsh.model.addPhysicalGroup(2, volumeTags, tag=1, name="Domain")

    if popup:
        gmsh.fltk.run()
    else:
        mesh = mhs.gmshModelToMesh( gmsh.model )
        gmsh.finalize()
        return mesh

if __name__=="__main__":
    getMeshPhysical(popup = True)
