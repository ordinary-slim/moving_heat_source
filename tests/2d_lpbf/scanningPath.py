import yaml
from numpy import round

params = {}
with open("input.yaml", 'r') as paramsFile:
    params = yaml.safe_load(paramsFile)

partLen = params["part"]
layerThickness = params["layerThickness"]
speed = max(params["HeatSourceSpeed"])

def writeGcode(gcodeFile="Path.gcode", nLayers=-1):
    gcodeLines = []
    # Start in -Y
    X = - partLen[0] / 2
    Y = 0.0
    Z = 0.0
    E = 0.0
    numLayers = round(partLen[1] / layerThickness).astype(int)
    if nLayers > 0:
        numLayers = min( numLayers, nLayers )
    
    gcodeLines.append( "G0 F{}".format(speed, X, Y, Z) )
    for ilayer in range(numLayers):
        Y = layerThickness * (ilayer + 1)
        E += 0.1
        gcodeLines.append( "G0 X{} Y{:.2f}".format(X, Y) )

        gcodeLines.append( "G4 P{}".format( params["interLayerDelay"] / 2 ) )
        gcodeLines.append( "G4 P{} R1".format( params["interLayerDelay"] / 2 ) )

        X = -X
        gcodeLines.append( "G1 X{} E{:.2f}".format(X, E) )

    with open(gcodeFile, 'w') as gfile:
        gfile.writelines( [line+"\n" for line in gcodeLines] )

if __name__=="__main__":
    writeGcode("Path.gcode")
