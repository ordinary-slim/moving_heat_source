import matplotlib.pyplot as plt
import numpy as np
import os

class myPlotHandler:
    def __init__(self, m, Tmin=-1, Tmax=+1, 
            pauseTime=0.25,
            figureFolder="",
            shortDescription=""):
        self.left = min(m.pos)
        self.right = max(m.pos)
        self.L = self.right - self.left
        self.currentTime = -1
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.range = abs(self.Tmax - self.Tmin)
        self.figCleared = False
        self.shortDescription = shortDescription
        self.figureFolder = figureFolder
        self.pauseTime = pauseTime
        self.mhsPlotted = False

    def clf(self, mesh):
        plt.clf()
        self.figCleared = True
        self.mhsPlotted = False
        #plt.plot(mesh.pos, np.zeros( len(mesh.pos) ), '-o')

    def plotProblem( self, p, label="solution",
            linestyle='-',
            updateLims=True,
            plotMhs=True,
            ):
        plt.plot( p.mesh.pos, p.solution, linestyle=linestyle, label=label );
        plt.xlim( self.left, self.right );
        #get max min temperatures
        Tmin = min( p.solution )
        Tmax = max( p.solution )

        if updateLims:
            self.Tmin  = min(self.Tmin, Tmin)
            self.Tmax  = max(self.Tmax, Tmax)
            self.range = max((self.Tmax-self.Tmin), self.range)
        if plotMhs and not(self.mhsPlotted):
            x0 = p.mhs.currentPosition[0]
            plt.plot(x0, self.Tmin, '-o',
                    color="red", label="$x_0$")
            plt.axvline(x=x0, linestyle='--', linewidth=0.5, color="red")
            self.mhsPlotted = True
        self.currentTime = p.time


        plt.ylim( self.Tmin - self.range*0.1, self.Tmax + self.range*0.2)

        plt.annotate("$t = {}$".format(str(round(p.time, 2))), (self.left+0.5*self.L, self.Tmin+0.05*self.range),
                );
        plt.xlabel(r"x $[mm]$")
        plt.ylabel(r"T $[{}^\circ C]$")
        yticks = list(np.linspace( self.Tmin, round(self.Tmax/1000, 1)*1000.0, 4 ))
        plt.yticks(yticks)
        plt.legend();

    def pause( self ):
        plt.pause( self.pauseTime )

    def save( self):
        if not self.figureFolder:
            self.figureFolder = "figures_" + self.shortDescription

        figureName = self.figureFolder + "/" + "t_{}.png".format( str(round(self.currentTime, 3)).replace(".", "_") )
        plt.savefig( figureName, dpi=300,
                bbox_inches='tight', )
