clear all;

inputdset = "inputWorkspace.mat";
myfrfscheme = FrfScheme(inputdset);
referenceSol = FrfScheme(inputdset);
myscheme = MyScheme(inputdset);

referenceSol.dt = 0.01;

myfrfscheme.preLoopAssembly();
referenceSol.preLoopAssembly();
myscheme.preLoopAssembly();


figure('Position', [100 100 900 600])
while myscheme.t < myscheme.Tfinal
    myfrfscheme.iterate();
    myscheme.iterate();
    while referenceSol.t < myscheme.t
        referenceSol.iterate();
    end
    % PLOT
    hold off
    plot(myfrfscheme.xpos, myfrfscheme.U, ...
        'DisplayName', "FRF", ...
        "LineWidth", 2);
    hold on
    plot(myscheme.pos+myscheme.t*myscheme.speed, myscheme.Upos, ...
        'DisplayName', "My scheme", ...
            "LineWidth", 2);
    plot(referenceSol.xpos, referenceSol.U, '--', ...
        'DisplayName', "Reference", ...
        "LineWidth", 1.5);
    xlim([myscheme.leftBound, myscheme.rightBound]);
    legend('Location', 'best')
    title(sprintf("Step size = %d R",myscheme.adimDt))
    set(gca, 'FontSize', 24)
    pause(0.4)
    % END PLOT
end