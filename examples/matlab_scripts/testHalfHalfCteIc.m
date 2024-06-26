clear all;

inputdset = "cteIc.mat";

S = load(inputdset, "leftBound", "rightBound", "power", ...
    "efficiency", "radius", "cutoffRadius", "x0", ...
    "speed", "rho", "cp", "k", "dt", "meshDensity", "Tfinal", "icX", "icXi");
S.icXi = @(xi, t) S.icX(xi+S.speed*t, t);
S.dt = 0.5;
S.meshDensity = 0.2;
S.k = 1.0;
S.cp = 1.0;
S.rho = 1.0;
S.power = 0.0;
S.Tfinal = 2.0;
S.dt = 1.0;
frfscheme = FrfScheme(S);
referenceSol = FrfScheme(S);
halfhalf = HalfHalfScheme(S);
S.meshDensity = 2*S.meshDensity;
halfhalfRef = HalfHalfScheme(S);

referenceSol.dt = 0.01;
referenceSol.meshDensity = 10;
referenceSol.initialize(S.icX);

frfscheme.preLoopAssembly();
referenceSol.preLoopAssembly();
halfhalf.preLoopAssembly();
halfhalfRef.preLoopAssembly();


figure('Position', [100 100 1200 900])
while halfhalf.t < halfhalf.Tfinal-1e-7
    frfscheme.iterate();
    halfhalf.iterate();
    halfhalfRef.iterate();
    while referenceSol.t < halfhalf.t
        referenceSol.iterate();
    end
    % PLOT
    hold off
    plot(frfscheme.xpos, frfscheme.U, ...
        'DisplayName', "FRF", ...
        "LineWidth", 2);
    hold on    
    plot(halfhalf.pos, halfhalf.Upos, ...
        'DisplayName', "Half half DD", ...
            "LineWidth", 2);
    plot(halfhalfRef.pos, halfhalfRef.Upos, ...
        'DisplayName', "Half half DD, refined in space", ...
            "LineWidth", 2);
    xline(halfhalf.xInterface, ...
        'DisplayName', "$\Gamma$")
    plot(referenceSol.xpos, referenceSol.U, '--', ...
        'DisplayName', "Reference", ...
        "LineWidth", 1.5);
    xlim([halfhalf.leftBound, halfhalf.rightBound]);
    % Support of heat source
    % myYlim = ylim();
    % rectangleColor = [1 0 0 0.05];
    % rectangle('Position', [frfscheme.x0-frfscheme.cutoffRadius myYlim(1) 2*frfscheme.cutoffRadius myYlim(2)-myYlim(1)], ...
    %     'FaceColor', [1 0 0 0.05], ...
    %     "LineStyle", "none")
    % scatter(NaN, NaN, 64, ...
    %     "Marker", 's', ...
    %     "MarkerEdgeColor", rectangleColor(1:3), "MarkerEdgeAlpha", rectangleColor(end), ...
    %     "MarkerFaceColor", rectangleColor(1:3), "MarkerFaceAlpha", rectangleColor(end), ...
    %     "DisplayName", "Heat source")
    delete(findall(gcf,'type','annotation'))
    dim = [.2 .5 .3 .3];
    timeString = sprintf("t = %.1f", halfhalf.t);
    annotation('textbox',dim,'String',timeString,'FitBoxToText','on', ...
        'Interpreter', 'latex', 'FontSize', 24);
    legend('Location', 'best', 'FontSize', 24, 'Interpreter', 'latex');
    title(sprintf("$\\Delta t$ = %.1f, h = %1.f",halfhalf.dt, halfhalf.h), ...
        'FontSize', 32, ...
        'Interpreter', 'latex')
    set(gca, 'FontSize', 24)
    set(gca, 'TickLabelInterpreter', 'latex')
    pause(0.4)
    % END PLOT
end
