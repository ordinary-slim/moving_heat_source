import MovingHeatSource as mhs
import numpy as np
from MovingHeatSource.adaptiveStepper import meshBox

nsteps=5

def deactivateBelowSurface(p, surfaceZ = 0):
    nels = p.domain.mesh.nels
    activeEls = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElement( ielem )
        if (e.getCentroid()[2] < surfaceZ):
            activeEls.append( ielem )
    substrateEls = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, activeEls )
    p.domain.setActivation( substrateEls )

def setAdimR( adimR, input ):
    r = input["radius"]
    HeatSourceSpeedX = max( abs(input["HeatSourceSpeedX"]), abs(input["advectionSpeedX"]))
    HeatSourceSpeedY = max( abs(input["HeatSourceSpeedY"]), abs(input["advectionSpeedY"]))
    HeatSourceSpeedZ = max( abs(input["HeatSourceSpeedZ"]), abs(input["advectionSpeedZ"]))
    speed  = np.linalg.norm( np.array( [HeatSourceSpeedX, HeatSourceSpeedY, HeatSourceSpeedZ] ) )
    return (adimR * r / speed)

def run():
    inputFile = "input.yaml"
    boxDomain = np.array([-0.005, 0.005, -0.025, +0.025, 0.0, 0.02])

    # read input
    problemInput = mhs.readInput( inputFile )
    mdwidth = 2.5e-3
    mdheight = 2.5e-3
    problemInput["heatSourceHeight"] = mdheight
    problemInput["heatSourceWidth"] = mdwidth

    fixedProblemInput = dict( problemInput )

    # Mesh
    meshFixed  = meshBox( boxDomain, elSize=0.01 )

    # mhs.Problem params
    # set dt
    dt = 1
    fixedProblemInput["dt"] = dt

    p           = mhs.Problem(meshFixed, fixedProblemInput, caseName="case")

    for p in [p]:
        deactivateBelowSurface( p, surfaceZ=0.01 )

    # Set up printer
    printerFRF = mhs.Printer( p, mdwidth, mdheight/2, mdheight/2 )

    for _ in range(nsteps):
        # Setup print
        p1 = p.mhs.position
        p2 = p1 + p.mhs.speed * dt
        # Print
        printerFRF.deposit( p1, p2, p.domain.activeElements )
        p.mhs.markHeatedElements(p1, p2)
        # Apply convection to all external boundaries
        p.setConvection()

        p.iterate()

        p.writepos()

def test():
    run()
    # COMPARISON
    allFilesSame = True
    for ifile in range(nsteps):
        newds = "post_case/case_{}.vtu".format( ifile+1 )
        refds = "post_case_reference/case_{}.vtu".format( ifile+1 )
        if not(mhs.meshio_comparison( newds, refds )):
               allFilesSame = False
               break
    assert allFilesSame

if __name__=="__main__":
    test()
