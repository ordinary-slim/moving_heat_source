clear all;

inputdset = "neumannLeft.mat";

S = load(inputdset, "leftBound", "rightBound", "power", ...
    "efficiency", "radius", "cutoffRadius", "x0", ...
    "speed", "rho", "cp", "k", "dt", "meshDensity", "Tfinal", "icX", "icXi", "neumannFluxLeft", "neumannFluxRight");
S.neumannFluxLeft = 50;
S.neumannFluxRight = 50.0;
S.dt = 0.2;
S.meshDensity = 2;
S.k = 1.0;
S.power = 0.0;
frfscheme = FrfScheme(S);
referenceSol = FrfScheme(S);
myscheme = MyScheme(S);
myschemeRef = MyScheme(S);
myschemeRefRef = MyScheme(S);

myschemeRef.dt = S.dt / 2;
myschemeRefRef.dt = S.dt / 4;
% myschemeRef.meshDensity = S.meshDensity * 2;
myschemeRef.minLengthSubdomainX = myscheme.getLengthSubdomainX();
myschemeRefRef.minLengthSubdomainX = myscheme.getLengthSubdomainX();
myschemeRef.initialize( S.icX, S.icXi );
myschemeRefRef.initialize( S.icX, S.icXi );

referenceSol.dt = 0.01;

frfscheme.preLoopAssembly();
referenceSol.preLoopAssembly();
myscheme.preLoopAssembly();
myschemeRef.preLoopAssembly();
myschemeRefRef.preLoopAssembly();

figure('Position', [100 100 1400 900])
while myscheme.t < myscheme.Tfinal-1e-7
    frfscheme.iterate();
    myscheme.iterate();
    while myschemeRef.t < myscheme.t
        myschemeRef.iterate();       
    end
    while myschemeRefRef.t < myscheme.t
        myschemeRefRef.iterate();       
    end
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
    plot(myschemeRef.pos+myschemeRef.t*myschemeRef.speed, myschemeRef.Upos, ...
        'DisplayName', "My scheme, refined in time", ...
            "LineWidth", 2);
    plot(myschemeRefRef.pos+myschemeRefRef.t*myschemeRefRef.speed, myschemeRefRef.Upos, ...
        'DisplayName', "My scheme, refined twice in time", ...
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
    title(sprintf("$\\Delta t$ = %.1f, h = %.1f",myscheme.dt, myscheme.h), ...
        'FontSize', 32, ...
        'Interpreter', 'latex')
    set(gca, 'FontSize', 24)
    set(gca, 'TickLabelInterpreter', 'latex')
    pause(0.4)
    % END PLOT
end