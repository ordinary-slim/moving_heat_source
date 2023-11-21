import meshzoo
import meshio
import numpy as np
from MovingHeatSource import readInput
import pdb

problemInput = readInput( "../input.yaml" )
power = problemInput["power"]
radius = problemInput["radius"]
rho = problemInput["material"]["rho"]
k = problemInput["material"]["conductivity"]
cp = problemInput["material"]["specific_heat"]
Tenv = problemInput["environmentTemperature"]
speed = np.linalg.norm( problemInput["HeatSourceSpeed"] )
kappa = k / rho / cp

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
    return points, cells

def solveNguyen(points, time):
    T = np.zeros( points.shape[0] )
    N = 6*np.sqrt(3) * power / (rho*cp*np.power(np.pi, 1.5))
    taoMax = time
    taoSamples = np.linspace( 0, taoMax, num=4000 )
    for inode, pos in enumerate(points):
        # NGuyen solution
        integrand = np.zeros( taoSamples.size )
        for idx, tao in enumerate(taoSamples):
            A1 = -3*((pos[0] + speed*(time-tao))**2 + pos[1]**2 + pos[2]**2)
            A1 /= (12*kappa*(time-tao) + radius**2)
            A1 = np.exp( A1 )
            integrand[idx] = A1 / np.power( (12*kappa*(time-tao) + radius**2), 1.5)
        integral = np.trapz( integrand, taoSamples )
        T[inode] = Tenv + N * integral
    return T

def main(time=0.01):
    elSize = radius/8
    # Mesh
    box = [-4*radius, +1*radius, 0.0, +2.0*radius, -elSize, 0.0]
    points, cells = meshBox( box, elSize )
    # Solve
    T = solveNguyen( points, time )
    # Posrun/3d_lpbft
    fileName = "nguyen.vtk"
    dataSet = meshio.Mesh(
            points,
            [("hexahedron", cells)],
            point_data={
                "T" : T,
                }
            )
    dataSet.write( fileName )

if __name__=="__main__":
    main()
