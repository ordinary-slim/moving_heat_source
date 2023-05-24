clear;
inputdset = "myDDPrototype.mat";
params = load(inputdset, "leftBound", "rightBound", "power", ...
    "efficiency", "radius", "cutoffRadius", "x0", ...
    "speed", "rho", "cp", "k", "dt", "meshDensity", "Tfinal", "icX", "icXi");

% PARAMS
params.Tenv = 25;
params.power = 100.0;
params.meshDensity =2;
params.Tfinal = 3;
% params.k = 0.0;
params.leftBound = -25;
params.rightBound = 25;
params.x0 = -15;
params.radiusSubdomain = 5*params.radius;
setDt( params, 2 );
params.ic = @(x) params.Tenv*ones(size(x));

bestSchemeEver = MyDDScheme( params );

tol = 1e-7;
figure
while params.Tfinal-tol > bestSchemeEver.problemPart.time
    bestSchemeEver.iterate();
    %% plot
    hold off
    plot(bestSchemeEver.problemPart.mesh.posFixed, bestSchemeEver.problemPart.U, 'r')
    hold on
    pause(0.25)
end

function [dt] = setDt( S, adimR)
    dt = S.radius / S.speed * adimR; 
end