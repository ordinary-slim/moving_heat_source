import yaml

params = {}
with open("input.yaml", 'r') as paramsFile:
    params = yaml.safe_load(paramsFile)

partLen = params["part"]
layerThickness = params["layerThickness"]
speed = max(params["HeatSourceSpeed"])
gcodeFile = params["path"]

def writeGcode(nLayers=None):
    gcodeLines = []
    # Start in -X
    X = - partLen[0] / 2
    Y = 0.0
    Z = 0.0
    E = 0.0
    numLayers = int(partLen[2] / layerThickness)
    if nLayers is not None:
        numLayers = min( numLayers, nLayers )
    
    gcodeLines.append( "G0 F{}".format(speed, X, Y, Z) )
    for ilayer in range(numLayers):
        Z = layerThickness * (ilayer)
        E += 0.1
        gcodeLines.append( "G0 X{} Y{} Z{:.2f}".format(X, Y, Z) )
        X = -X
        gcodeLines.append( "G1 X{} E{:.2f}".format(X, E) )

    with open(gcodeFile, 'w') as gfile:
        gfile.writelines( [line+"\n" for line in gcodeLines] )

if __name__=="__main__":
    writeGcode()
