import MovingHeatSource as mhs
from MovingHeatSource.adaptiveStepper import AdaptiveStepper, TrackType
import numpy as np
import meshio

class CustomStepper(AdaptiveStepper):
    def shapeSubdomain( self ):
        '''
        At t^n, do things
        '''
        if not(self.nextTrack.type == TrackType.printing):
            return
        # OBB
        radius = self.pFixed.mhs.radius

        # compute front and sides
        sideRadius = self.adimMinRadius * radius
        adimBackRadius = min( self.adimMaxSubdomainSize, self.adimSubdomainSize )
        backRadius = max( adimBackRadius, self.adimMinRadius ) * radius
        zRadius    = 1
        xAxis      = self.nextTrack.getSpeed() / self.nextTrack.speed

        backRadiusObb = max(backRadius - radius, 0.0)
        p0 = self.pMoving.mhs.position - backRadiusObb*xAxis
        obb = mhs.MyOBB( p0, self.pMoving.mhs.position, 2*sideRadius, 2*zRadius )
        subdomainEls = self.pMoving.domain.mesh.findCollidingElements( obb )
        collidingElsBackSphere = self.pMoving.domain.mesh.findCollidingElements( p0, self.adimMinRadius*radius )
        collidingElsFrontSphere = self.pMoving.domain.mesh.findCollidingElements( self.pMoving.mhs.position, self.adimMinRadius*radius )
        subdomainEls += collidingElsBackSphere
        subdomainEls += collidingElsFrontSphere
        subdomain = mhs.MeshTag( self.pMoving.domain.mesh, self.pMoving.domain.mesh.dim, subdomainEls )
        self.pMoving.domain.intersect( subdomain )

    def writepos( self ):
        activeInExternal = self.pFixed.getActiveInExternal( self.pMoving, 1e-7 )
        self.pFixed.writepos(
            nodeMeshTags={
                "gammaNodes":self.pFixed.gammaNodes,
                "forcedDofs":self.pFixed.forcedDofs,
                "activeInExternal":activeInExternal,
                },
            cellMeshTags={
                "physicalDomain":self.physicalDomain,
                },
                          )
        self.pMoving.writepos(
            nodeMeshTags={ "gammaNodes":self.pMoving.gammaNodes, },
            )

def deactivateBelowSurface(p,
                           surfacey = 0):
    nels = p.domain.mesh.nels
    activeels = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElementGeometry( ielem )
        if (e.getCentroid()[1] < surfacey):
            activeels.append( ielem )
    substrateels = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, activeels )
    p.domain.setActivation( substrateels )

def readMesh(gmshFile):
    m = meshio.read( gmshFile )
    mDict = {}
    mDict["points"] = m.points

    cells = np.vstack( [cell.data for cell in m.cells] )
    mDict["cells"] = cells
    mDict["cell_type"] = "quad4"
    mDict["dimension"] = 2

    return mhs.Mesh( mDict )

def main():
    inputFile = "input.yaml"
    problemInput = mhs.readInput( inputFile )

    # read input
    fixedProblemInput = dict( problemInput )

    # Mesh
    meshFixed = readMesh("untitled.msh")

    pFixed   = mhs.Problem(meshFixed, fixedProblemInput, caseName="fixed")

    # Set path, deactivate
    deactivateBelowSurface( pFixed )

    myDriver = CustomStepper( pFixed,
                              maxAdimtDt=2,
                              threshold=0.025,
                               )

    #while not(pFixed.mhs.path.isOver(pFixed.time)):
    for _ in range(4):
        myDriver.iterate()

def test():
    main()
    reference = "post_fixed_reference/fixed_4.vtu"
    new = "post_fixed/fixed_4.vtu"
    assert mhs.meshio_comparison(reference, new)

if __name__=="__main__":
    test()
