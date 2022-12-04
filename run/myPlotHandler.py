import matplotlib.pyplot as plt
import numpy as np

class myPlotHandler:
    def __init__(self, m, Tmin=-1, Tmax=+1):
        self.left = min(m.pos)
        self.right = max(m.pos)
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.range = abs(self.Tmax - self.Tmin)
        self.figCleared = False

    def clf(self, mesh):
        plt.clf()
        self.figCleared = True
        #plt.plot(mesh.pos, np.zeros( len(mesh.pos) ), '-o')

    def plotProblem( self, p, label="solution",
            updateLims=True,
            plotMhs=True,
            ):
        plt.plot( p.mesh.pos, p.solution, label=label );
        plt.xlim( self.left, self.right );
        #get max min temperatures
        Tmin = min( p.solution )
        Tmax = max( p.solution )

        if updateLims:
            self.Tmin  = min(self.Tmin, Tmin)
            self.Tmax  = max(self.Tmax, Tmax)
            self.range = max((self.Tmax-self.Tmin), self.range)
        if plotMhs:
            plt.plot(p.mhs.currentPosition[0], self.Tmin, '-o',
                    color="red")


        plt.ylim( self.Tmin - self.range*0.1, self.Tmax + self.range*0.1)

        plt.annotate("time = " + str(round(p.time, 2)), (0.1, 0.9),
                xycoords='figure fraction');
        plt.legend();

