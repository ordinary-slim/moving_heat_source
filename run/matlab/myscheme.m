clc; clear all;
%% PARAMS
L = 50;
leftBound = -L/2;
rightBound = leftBound + L;

% HEAT SOURCE
power = 100.0;
efficiency = 1.0;
radius = 1;
xi0 = -5.0;
speed = 10;
% IC
icXi = @(xi, t) 25*ones(size(xi));
icX = @(x, t) 25*ones(size(x));
% power = 0.0;
% icX = @(x, t) max(0.0, 20-abs(x));
% icXi = @(xi, t) icX(xi+speed*t, t);
% power = 0.0;
% icX = @(x, t) sin(2*(x)/L*pi);
% icXi = @(xi, t) icX(xi+speed*t, t);

% MATERIAL
rho = 1;
cp = 1;
k = 1;


%% DISCRETIZATION
% TIME
tol = 1e-7;
t = 0.0;
dt = 1;
Tfinal = 2.0;
% SPACE
meshDen = 5;
h = 1/meshDen;
xi_nodesFun = @(t) (-speed*t+leftBound):h:(-speed*(t+dt) + rightBound);
x_nodesFun = @(t) (leftBound):h:(leftBound+speed*dt);
xipos = xi_nodesFun(t);
xpos  = x_nodesFun(t);
nnodesXi = size(xipos,2);
nnodesX  = size(xpos, 2);
nelsXi = nnodesXi-1;
nelsX = nnodesX-1;
xi_connectivity = [(1:nelsXi)', (2:(nelsXi+1))'];
x_connectivity = [(1:nelsX)', (2:(nelsX+1))'];

nnodes = nnodesXi + nnodesX;
nels   = nelsXi + nelsX;

% Build global numbering
%numberingSubproblem(i) gives global numbering of node i of Subproblem
numberingX  = 1:nnodesX;
numberingXi = (nnodesX+1):(nnodesX+nnodesXi);

%% Initializing solution
Uxi = icXi(xipos, 0)';
Ux = icX(xpos, 0)';
U = [Ux; Uxi];
massX = sparse(nnodesX, nnodesX);
massXi = sparse(nnodesXi, nnodesXi);
diffusionX = sparse(nnodesX, nnodesX);
diffusionXi = sparse(nnodesXi, nnodesXi);
advectionXi = sparse(nnodesXi, nnodesXi);
pulse = zeros([nnodesXi, 1]);

% FIXED SUBPROBLEM
% Assemble mass matrix
for iel=1:nelsX
    inodes = x_connectivity(iel, :);
    xposloc = xpos(inodes);
    h = xposloc(2) - xposloc(1);
    Mloc = [h/3, h/6; h/6, h/3];
    Kloc = [1/h, -1/h; -1/h, 1/h];
    massX(inodes, inodes) = ...
        massX(inodes, inodes) + Mloc;
    diffusionX(inodes, inodes) = ...
        diffusionX(inodes, inodes) + Kloc;
end
% XI SUBPROBLEM
% Assemble advection and mass mat
for iel=1:nelsXi
    inodes = xi_connectivity(iel, :);
    xiposloc = xipos(inodes);
    h = xiposloc(2) - xiposloc(1);
    Mloc = [h/3, h/6; h/6, h/3];
    Aloc = [-1/2, 1/2; -1/2, 1/2];
    Kloc = [1/h, -1/h; -1/h, 1/h];
    massXi(inodes, inodes) = ...
        massXi(inodes, inodes) + Mloc;
    advectionXi(inodes, inodes) = ...
        advectionXi(inodes, inodes) + Aloc;
    diffusionXi(inodes, inodes) = ...
        diffusionXi(inodes, inodes) + Kloc;
end

%% Time loop
iter = 0;
figure
while t < Tfinal
    iter = iter + 1;
    "iter " + iter
    xipos = xi_nodesFun(t);
    xpos = x_nodesFun(t);
    t = t + dt;
    %% Assembly    
    lhs = sparse(nnodes, nnodes);
    rhs = zeros([nnodes, 1]);
    pulse = zeros([nnodesXi, 1]);
    % Assemble heat source into rhs xi
    for iel=1:nelsXi
        inodes = xi_connectivity(iel, :);
        xiposloc = xipos(inodes);
        h = xiposloc(2) - xiposloc(1);
        rloc = h*powerDensity(xiposloc, xi0, power, efficiency, radius);
        pulse(inodes) = pulse(inodes) + rloc';
    end
    % Assemble subproblems into system
    rhs(numberingX) = rhs(numberingX) + rho*cp*( massX*Ux / dt );
    rhs(numberingXi) = rhs(numberingXi) + rho*cp*( massXi*Uxi / dt ) + pulse;
    lhs(numberingX, numberingX) = lhs(numberingX, numberingX) + rho*cp*(massX/dt) + k*diffusionX;
    lhs(numberingXi, numberingXi) = lhs(numberingXi, numberingXi) + rho*cp*(massXi/dt - speed*advectionXi) + k*diffusionXi;

    % Assemble Neumann condition interface
    % inodeBounX = x_connectivity( nelsX, 2);
    % inodesXi = xi_connectivity( 1, :);
    % xiposloc = xipos(inodesXi);
    % h = xiposloc(2) - xiposloc(1);
    % Kboun = [-1/h, 1/h];
    % lhs(numberingX(inodeBounX), numberingXi(inodesXi)) = ...
    %     lhs(numberingX(inodeBounX), numberingXi(inodesXi)) - k*Kboun;
    % % TODO, fix here

    % Assemble Dirichlet interface RIGHT
    % inodeDirichletXi = xi_connectivity( 1, 1);
    % inodeDirichletX = x_connectivity( nelsX, 2 );
    % lhs(numberingXi(inodeDirichletXi), :) = 0.0;
    % lhs(numberingXi(inodeDirichletXi), numberingX(inodeDirichletX)) = -1.0;
    % lhs(numberingXi(inodeDirichletXi), numberingXi(inodeDirichletXi)) = 1.0;
    % rhs(numberingXi(inodeDirichletXi)) = 0.0;

    % Assemble Neumann condition interface RIGHT
    inodeBounXi_Neumann = xi_connectivity( 1, 1);
    inodesX_Neumann = x_connectivity( nelsX, :);
    xposloc = xpos(inodesX_Neumann);
    inodeBounXi_Neumann_global = numberingXi( inodeBounXi_Neumann );
    inodesX_Neumann_global = numberingX( inodesX_Neumann );
    h = xposloc(2) - xposloc(1);
    Kboun = -[-1/h, 1/h];
    lhs(inodeBounXi_Neumann_global, inodesX_Neumann_global) = ...
        lhs(inodeBounXi_Neumann_global, inodesX_Neumann_global) - k*Kboun;
    % TODO, fix here
    % Assemble Dirichlet interface LEFT
    inodeDirichletXi_global = numberingXi( xi_connectivity( 1, 1) );
    inodeDirichletX_global = numberingX( x_connectivity( nelsX, 2 ) );
    lhs(inodeDirichletX_global, :) = 0.0;
    lhs(inodeDirichletX_global, inodeDirichletX_global) = 1.0;
    lhs(inodeDirichletX_global, inodeDirichletXi_global) = -1.0;
    rhs(inodeDirichletX_global) = 0.0;
    %% Solve
    U = lhs \ rhs;

    %% pos-iteration
    posSol = [xpos(1:end-1)-speed*t, xipos];
    sol = [U(numberingX(1:end-1)); U(numberingXi)];
    Ux = U(numberingX);
    Uxi = interp1(posSol, sol, xi_nodesFun(t))';%gotta work on this
    % Uxi = exactSolXi(xi_nodesFun(t), t)';
    
    minY = min(U);
    maxY = max(U);
    range = maxY - minY;
    hold off
    plot(posSol+speed*t, sol, 'DisplayName', "myScheme");
    legend()
    % xlim([leftBound-speed*(Tfinal+3*dt), rightBound+speed*dt]);
    xlim([leftBound, rightBound]);
    ylim([minY-0.1*range, maxY+0.1*range]);
    pause(0.4);
end

%% HEAT SOURCE
function [pd] = powerDensity(xi, xi0, power, efficiency, radius)
    if (abs(xi-xi0)>3*radius)
        pd = 0.0;
    else
        pd = 2*(power*efficiency) / pi / radius^2 * exp( - 2*(xi - xi0).^2/radius^2);
    end
end