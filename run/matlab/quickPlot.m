clear all;

inputdset = "gaussianTest.mat";

S = load(inputdset, "leftBound", "rightBound", "power", ...
    "efficiency", "radius", "cutoffRadius", "x0", ...
    "speed", "rho", "cp", "k", "dt", "meshDensity", "Tfinal", "icX", "icXi");
S.x0 = 0.0;
S.power = 100;
S.radius = 2;
frfscheme = FrfScheme(S);

R = S.radius;

xrange = linspace( -4*R, 4*R );
pd = frfscheme.powerDensity( xrange, 0.0 );

figure('Position', [100 100 800 600])
plot( xrange, pd, "Color", "red", 'LineWidth', 2)
range = max(pd);
ylim( [0, 1.1*range] )
xlim( [xrange(1), xrange(end)] )
xlabel("x")
ylabel("Power density")
myXticks = xticks();
% myXticks = sort( [myXticks, -2*S.radius, 2*S.radius] );
xticks(myXticks)
title(sprintf("R = %g,   power = %g", frfscheme.radius, frfscheme.power ))
set(gca, 'FontSize', 24)
grid on