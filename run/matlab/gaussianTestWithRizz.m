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

adimDt = 0.5;
adimDomainSize = 12;
adimPad = 0.1*adimDomainSize;
params.x0 = -20;
params.radiusSubdomain = adimDomainSize*params.radius;
params.dt = setDt( params, adimDt );
params.pad = adimPad*params.radius;

fineparams = params;
fineparams.meshDensity = 4;
fineparams.dt = setDt( params, 0.2 );

bestSchemeEver = MyDDScheme( params );
bestSchemeEver.label = "My scheme";
schwarz = SchwarzDDScheme( params );
frfscheme = FRFScheme( params );
myoldscheme = MyScheme( params );
finefrfscheme = FRFScheme( fineparams );

tol = 1e-7;
figure('Position', [200 100 1200 900])

while params.Tfinal-tol > bestSchemeEver.problemPart.time
    % if (bestSchemeEver.problemPart.time < 0.4 - tol)
    %     adimDt = 0.5;
    %     adimDomainSize = 1;
    if (bestSchemeEver.problemPart.time < 1 - tol)
        adimDt = 1;
        adimDomainSize = 2;
    elseif (bestSchemeEver.problemPart.time < 1.5 - tol)
        adimDt = 2.5;
        adimDomainSize = 4;
    elseif (bestSchemeEver.problemPart.time < 3 - tol)
        adimDt = 5;
        adimDomainSize = 7;
    else
        adimDt = 10;
        adimDomainSize = 12;
    end
    bestSchemeEver.setDt( params.radius / params.speed * adimDt );
    bestSchemeEver.setRadiusSubdomain( adimDomainSize * params.radius);
    
    bestSchemeEver.iterate();
    % schwarz.iterate();

    % myoldscheme.iterate();
    % while finefrfscheme.getTime < bestSchemeEver.getTime()-tol
    %     finefrfscheme.iterate();
    % end
    while frfscheme.getTime < bestSchemeEver.getTime-tol
        frfscheme.iterate();
    end
    %% plot
    hold off
    plot(bestSchemeEver.problemPart.mesh.posFixed, bestSchemeEver.problemPart.U, ...
         "DisplayName", bestSchemeEver.label + ", variable $\Delta t$", "LineWidth", 2)
    hold on
    % plot(schwarz.problemPart.mesh.posFixed, schwarz.problemPart.U, ...
        % "DisplayName", "Schwarz", "LineWidth", 2)
    
    % plot(schwarz.problemHS.mesh.posFixed, schwarz.problemHS.U, ...
         % "DisplayName", "Schwarzie mv", "LineWidth", 2)
    plot(frfscheme.mesh.posFixed, frfscheme.problem.U, ...
         "DisplayName", "FRF scheme, $\Delta t = 0.5 \mathcal{R}$", "LineWidth", 2)
    % plot(finefrfscheme.mesh.posFixed, finefrfscheme.problem.U, '--', ...
    %     'DisplayName', "Reference", ...
    %     "LineWidth", 1.5);
    % plot(myoldscheme.pos+myoldscheme.t*myoldscheme.speed, myoldscheme.Upos, ...
    %     'DisplayName', "Almost everywhere MRF", ...
    %         "LineWidth", 2);
    
    % PLOT INTERFACE
    % schwarz.plotInterface()
    bestSchemeEver.plotInterface()
    
    delete(findall(gcf,'type','annotation'))
    dim = [.45 0.0 .3 .3];
    timeString = sprintf("t = %.1f", bestSchemeEver.getTime());
    annotation('textbox',dim,'String',timeString,'FitBoxToText','on', ...
        'Interpreter', 'latex', 'FontSize', 24);
    xlim([params.leftBound, params.rightBound]);
    legend('Location', 'best', 'FontSize', 24, 'Interpreter', 'latex');
    title(sprintf("Current $\\Delta t$ = %g = %g $\\mathcal{R}$, h = %g, subdomain = %g $\\mathcal{R}$, V = %g", ...
        bestSchemeEver.problemPart.dt, adimDt, 1/params.meshDensity, adimDomainSize, params.speed), ...
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
