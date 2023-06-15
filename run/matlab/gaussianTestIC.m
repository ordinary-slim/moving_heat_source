clear;
inputdset = "myDDPrototype.mat";
params = load(inputdset, "leftBound", "rightBound", "power", ...
    "efficiency", "radius", "cutoffRadius", "x0", ...
    "speed", "rho", "cp", "k", "dt", "meshDensity", "Tfinal", "icX", "icXi");

% PARAMS
params.Tenv = 25;
params.power = 100.0;
params.meshDensity = 2;
params.Tfinal = 4;
params.leftBound = -25;
params.rightBound = +25;
params.ic = @(x) params.Tenv*ones(size(x));

adimDt = 10;
adimDomainSize = 11;
adimPad = 0.1*adimDomainSize;
params.x0 = -10;
params.radiusSubdomain = adimDomainSize*params.radius;
params.dt = setDt( params, adimDt );
params.pad = adimPad*params.radius;

fineparams = params;
fineparams.meshDensity = 2;
fineparams.dt = setDt( params, 0.2 );

bestSchemeEver = MyDDScheme( params );
frfscheme = FRFScheme( params );
myoldscheme = MyScheme( params );
finefrfscheme = FRFScheme( fineparams );

tol = 1e-7;
figure('Position', [200 100 1200 900])

while params.Tfinal-tol > bestSchemeEver.problemPart.time
    bestSchemeEver.iterate();
    % frfscheme.iterate();
    % myoldscheme.iterate();
    % while finefrfscheme.getTime < bestSchemeEver.getTime()
    %     finefrfscheme.iterate();
    % end
    %% plot
    hold off
    plot(bestSchemeEver.problemPart.mesh.posFixed, bestSchemeEver.problemPart.U, ...
         "DisplayName", "$U^{n+1}$", "LineWidth", 2)
    hold on
    % scatter(bestSchemeEver.problemPart.mesh.posFixed(logical(bestSchemeEver.problemPart.mesh.activeNodes)), ...
    %     bestSchemeEver.problemPart.Uprev, ...
    %      24, "DisplayName", "$U^{n}$ fixed subdomain", "LineWidth", 2)
    % scatter(bestSchemeEver.problemHS.mesh.posFixed(logical(bestSchemeEver.problemHS.mesh.activeNodes)), ...
    %     bestSchemeEver.problemHS.Uprev, ...
    %      24, "DisplayName", "$U^{n}$ moving subdomain", "LineWidth", 2)

    bestSchemeEver.plotInterface();

    delete(findall(gcf,'type','annotation'))
    dim = [.45 0.0 .3 .3];
    timeString = sprintf("t = %.1f", bestSchemeEver.getTime());
    annotation('textbox',dim,'String',timeString,'FitBoxToText','on', ...
        'Interpreter', 'latex', 'FontSize', 24);
    xlim([params.leftBound, params.rightBound]);
    legend('Location', 'best', 'FontSize', 24, 'Interpreter', 'latex');
    title(sprintf("$\\Delta t$ = %g = %g $\\mathcal{R}$, h = %g, subdomain = %g $\\mathcal{R}$, V = %g", ...
        bestSchemeEver.problemPart.dt, adimDt, 1/params.meshDensity, ...
        bestSchemeEver.currRadiusSubdomain/params.radius, params.speed), ...
        'FontSize', 32, ...
        'Interpreter', 'latex')
    set(gca, 'FontSize', 24)
    set(gca, 'TickLabelInterpreter', 'latex')
    xlabel("x", "Interpreter", "latex")
    ylabel("u", "Interpreter", "latex")
    grid on
    pause(0.25)
end

function [dt] = setDt( S, adimR)
    dt = S.radius / S.speed * adimR;
end
