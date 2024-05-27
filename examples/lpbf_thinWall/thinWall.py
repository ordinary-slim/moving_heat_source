import gmsh
import numpy as np
import MovingHeatSource as mhs
import yaml
import sys

params = {}
with open("input.yaml", 'r') as paramsFile:
    params = yaml.safe_load(paramsFile)
partLen = params["partLen"]
substrateLen = params["substrateLen"]
layerThickness = params["layerThickness"]
radiusHeatSource = params["radius"]
fineElementSize = min( layerThickness, radiusHeatSource) / 2
bLayerEls = 4

def basicProgression( numEls ):
    numElements = [1]*numEls
    heights = np.cumsum( np.linspace(0.5, 1, numEls) )
    heights /= heights[-1]
    return numElements, heights

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

def getGmshModel(popup=False):
    # Bot surface part
    gmsh.model.geo.addPoint(-partLen[0]/2, -partLen[1]/2, 0, 0.1, 1)
    gmsh.model.geo.addPoint(+partLen[0]/2, -partLen[1]/2, 0, 0.1, 2)
    gmsh.model.geo.addPoint(+partLen[0]/2, +partLen[1]/2, 0, 0.1, 3)
    gmsh.model.geo.addPoint(-partLen[0]/2, +partLen[1]/2, 0, 0.1, 4)
    #
    smallLine1 = gmsh.model.geo.addLine(1, 2, 1)
    bigLine1 = gmsh.model.geo.addLine(2, 3, 2)
    smallLine2 = gmsh.model.geo.addLine(3, 4, 3)
    bigLine2 = gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([bigLine2, smallLine1, bigLine1, smallLine2], 1)
    tagBotSurfacePart = 1
    gmsh.model.geo.addPlaneSurface([1], tagBotSurfacePart)
    gmsh.model.geo.mesh.setRecombine(2,tagBotSurfacePart) 
    # One element per layer
    nelsPartZ = partLen[2] / fineElementSize
    topExtrusion = gmsh.model.geo.extrude([(2, tagBotSurfacePart)], 0, 0, +partLen[2], numElements=[nelsPartZ], recombine=True)
    #
    nelsSubstrateZ = int(substrateLen[2] / layerThickness / 1.5)
    numElements, heights = boundaryLayerProgression( substrateLen[2], bLayerEls, fineElementSize, coarseElFactor=4 )
    botExtrusion = gmsh.model.geo.extrude([(2, tagBotSurfacePart)], 0, 0, -substrateLen[2], numElements=numElements, heights=heights, recombine=True)
    #
    nelsXLinePart = int( partLen[0] / fineElementSize )
    for line in [smallLine1, smallLine2]:
        gmsh.model.geo.mesh.setTransfiniteCurve(line, nelsXLinePart+1)
    nelsYLinePart = int(partLen[1] / fineElementSize)
    for line in [bigLine1, bigLine2]:
        gmsh.model.geo.mesh.setTransfiniteCurve(line, nelsYLinePart+1)
    gmsh.model.geo.mesh.setTransfiniteSurface(tagBotSurfacePart)

    centralLongBotSurfaces = [botExtrusion[2], botExtrusion[4]]

    # Substrate X-extrusions
    extrusionLen = (substrateLen[0]-partLen[0])/2
    nelsXSubstrate = int(extrusionLen / fineElementSize)

    numElements, heights = boundaryLayerProgression( extrusionLen, bLayerEls, fineElementSize, coarseElFactor=4 )
    botSideExtrusion1 = gmsh.model.geo.extrude([centralLongBotSurfaces[0]], -extrusionLen, 0, 0, numElements=numElements, heights=heights, recombine=True)
    botSideExtrusion2 = gmsh.model.geo.extrude([centralLongBotSurfaces[1]], +extrusionLen, 0, 0, numElements=numElements, heights=heights, recombine=True)

    # Substrate Y-extrusions
    extrusionLen = (substrateLen[1]-partLen[1])/2
    numElements, heights = boundaryLayerProgression( extrusionLen, bLayerEls, fineElementSize, coarseElFactor=4 )
    botFrontExtrusionSurfaceTags = [botSideExtrusion2[3], botExtrusion[5], botSideExtrusion1[5]]
    botBackExtrusionSurfaceTags = [botSideExtrusion1[3], botExtrusion[3], botSideExtrusion2[5]]
    botFrontExtrusion = gmsh.model.geo.extrude(botFrontExtrusionSurfaceTags, +0.0, +extrusionLen, 0.0, numElements=numElements, heights=heights, recombine=True)
    botBackExtrusion = gmsh.model.geo.extrude(botBackExtrusionSurfaceTags, +0.0, -extrusionLen, 0.0, numElements=numElements, heights=heights, recombine=True)

    # Add physical group for domain
    volumeTags = []
    gmsh.model.geo.synchronize()
    for extrusion in [topExtrusion, botExtrusion, botFrontExtrusion, botBackExtrusion, botSideExtrusion1, botSideExtrusion2]:
        for dim, tag in extrusion:
            if dim==3:
                volumeTags.append( tag )

    gmsh.model.addPhysicalGroup(3, volumeTags, tag=1, name="Domain")

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    # Launch the GUI to see the results:
    if popup:
        gmsh.fltk.run()

    return gmsh.model

def getGmshModelUniform(popup=False):
    # Bot surface part
    gmsh.model.geo.addPoint(-partLen[0]/2, -partLen[1]/2, 0, 0.1, 1)
    gmsh.model.geo.addPoint(+partLen[0]/2, -partLen[1]/2, 0, 0.1, 2)
    gmsh.model.geo.addPoint(+partLen[0]/2, +partLen[1]/2, 0, 0.1, 3)
    gmsh.model.geo.addPoint(-partLen[0]/2, +partLen[1]/2, 0, 0.1, 4)
    #
    smallLine1 = gmsh.model.geo.addLine(1, 2, 1)
    bigLine1 = gmsh.model.geo.addLine(2, 3, 2)
    smallLine2 = gmsh.model.geo.addLine(3, 4, 3)
    bigLine2 = gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([bigLine2, smallLine1, bigLine1, smallLine2], 1)
    tagBotSurfacePart = 1
    gmsh.model.geo.addPlaneSurface([1], tagBotSurfacePart)
    gmsh.model.geo.mesh.setRecombine(2,tagBotSurfacePart) 
    # One element per layer
    nelsPartZ = partLen[2] / fineElementSize
    topExtrusion = gmsh.model.geo.extrude([(2, tagBotSurfacePart)], 0, 0, +partLen[2], numElements=[nelsPartZ], recombine=True)
    #
    nelsSubstrateZ = int(substrateLen[2] / fineElementSize)
    botExtrusion = gmsh.model.geo.extrude([(2, tagBotSurfacePart)], 0, 0, -substrateLen[2], numElements=[nelsSubstrateZ], recombine=True)
    #
    nelsXLinePart = int( partLen[0] / fineElementSize )
    for line in [smallLine1, smallLine2]:
        gmsh.model.geo.mesh.setTransfiniteCurve(line, nelsXLinePart+1)
    nelsYLinePart = int(partLen[1] / fineElementSize)
    for line in [bigLine1, bigLine2]:
        gmsh.model.geo.mesh.setTransfiniteCurve(line, nelsYLinePart+1)
    gmsh.model.geo.mesh.setTransfiniteSurface(tagBotSurfacePart)

    centralLongBotSurfaces = [botExtrusion[2], botExtrusion[4]]

    # Substrate X-extrusions
    extrusionLen = (substrateLen[0]-partLen[0])/2
    nelsXSubstrate = int(extrusionLen / fineElementSize)
    botSideExtrusion1 = gmsh.model.geo.extrude([centralLongBotSurfaces[0]], -extrusionLen, 0, 0, numElements=[nelsXSubstrate], recombine=True)
    botSideExtrusion2 = gmsh.model.geo.extrude([centralLongBotSurfaces[1]], +extrusionLen, 0, 0, numElements=[nelsXSubstrate], recombine=True)

    # Substrate Y-extrusions
    extrusionLen = (substrateLen[1]-partLen[1])/2
    numElsExtrusionYSubstrate = extrusionLen / fineElementSize
    botFrontExtrusionSurfaceTags = [botSideExtrusion2[3], botExtrusion[5], botSideExtrusion1[5]]
    botBackExtrusionSurfaceTags = [botSideExtrusion1[3], botExtrusion[3], botSideExtrusion2[5]]

    botFrontExtrusion = gmsh.model.geo.extrude(botFrontExtrusionSurfaceTags, +0.0, +extrusionLen, 0.0, numElements=[numElsExtrusionYSubstrate], recombine=True)
    botBackExtrusion = gmsh.model.geo.extrude(botBackExtrusionSurfaceTags, +0.0, -extrusionLen, 0.0, numElements=[numElsExtrusionYSubstrate], recombine=True)

    # Add physical group for domain
    volumeTags = []
    gmsh.model.geo.synchronize()
    for extrusion in [topExtrusion, botExtrusion, botFrontExtrusion, botBackExtrusion, botSideExtrusion1, botSideExtrusion2]:
        for dim, tag in extrusion:
            if dim==3:
                volumeTags.append( tag )

    gmsh.model.addPhysicalGroup(3, volumeTags, tag=1, name="Domain")

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    # Launch the GUI to see the results:
    if popup:
        gmsh.fltk.run()

    return gmsh.model


def generateMesh():
    gmsh.initialize()
    model = getGmshModel()
    mesh = mhs.gmshModelToMesh( model )
    gmsh.finalize()
    return mesh

if __name__=="__main__":
    gmsh.initialize()
    if "--uniform" in sys.argv:
        _ = getGmshModelUniform(popup=True)
    else:
        _ = getGmshModel(popup=True)
    gmsh.finalize()
