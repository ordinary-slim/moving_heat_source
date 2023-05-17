clear all;

inputdset = "neumannLeft.mat";

S = load(inputdset, "leftBound", "rightBound", "power", ...
    "efficiency", "radius", "cutoffRadius", "x0", ...
    "speed", "rho", "cp", "k", "dt", "meshDensity", "Tfinal", "icX", "icXi", "neumannFluxLeft", "neumannFluxRight");
S.neumannFluxLeft = 10;
S.neumannFluxRight = 0.0;
S.dt = 0.1;
S.meshDensity = 5;
S.k = 1.0;
S.power = 0.0;
frfscheme = FrfScheme(S);
referenceSol = FrfScheme(S);
myscheme = MyScheme(S);

referenceSol.dt = 0.01;

frfscheme.preLoopAssembly();
referenceSol.preLoopAssembly();
myscheme.preLoopAssembly();


figure('Position', [100 100 1200 900])
while myscheme.t < myscheme.Tfinal-1e-7
    frfscheme.iterate();
    myscheme.iterate();
    while referenceSol.t < myscheme.t
        referenceSol.iterate();
    end
    % PLOT
    hold off
    plot(frfscheme.xpos, frfscheme.U, ...
        'DisplayName', "FRF", ...
        "LineWidth", 2);
    hold on
    plot(myscheme.pos+myscheme.t*myscheme.speed, myscheme.Upos, ...
        'DisplayName', "My scheme", ...
            "LineWidth", 2);
    xline(myscheme.xInterface, ...
        'DisplayName', "My scheme, $\Gamma$")
    plot(referenceSol.xpos, referenceSol.U, '--', ...
        'DisplayName', "Reference", ...
        "LineWidth", 1.5);
    xlim([myscheme.leftBound, myscheme.rightBound]);
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
    timeString = sprintf("t = %.2f", myscheme.t);
    annotation('textbox',dim,'String',timeString,'FitBoxToText','on', ...
        'Interpreter', 'latex', 'FontSize', 24);
    legend('Location', 'best', 'FontSize', 24, 'Interpreter', 'latex');
    title(sprintf("$\\Delta t$ = %.2f",myscheme.dt), ...
        'FontSize', 32, ...
        'Interpreter', 'latex')
    set(gca, 'FontSize', 24)
    set(gca, 'TickLabelInterpreter', 'latex')
    %pause(0.4)
    % END PLOT
end