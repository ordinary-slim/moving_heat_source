clc; clear all;
%% PARAMS
L = 50;
leftBound = -L/2;
rightBound = leftBound + L;

% HEAT SOURCE
power = 100;
efficiency = 1.0;
radius = 1;
x0 = 0.0;
speed = 10;
% IC
icXi = @(xi, t) 25*ones(size(xi));
icX = @(x, t) 25*ones(size(x));

% MATERIAL
rho = 1;
cp = 1;
k = 1;


%% DISCRETIZATION
nelsXi = 80;%# elements
nnodesXi = nelsXi + 1;
nelsX = 20;%# elements
nnodesX = nelsX + 1;
nnodes = nnodesXi + nnodesX;
nels   = nelsXi + nelsX;
dx = L/nels;
% TIME
tol = 1e-7;
t = 0.0;
dt = 0.2;
Tfinal = 2.0;

%% MESHING
xi_nodesFun = @(t) linspace(-speed*t + leftBound,-speed*(t+dt) + rightBound, nnodesXi);
xi_connectivity = [(1:nelsXi)', (2:(nelsXi+1))'];
x_nodesFun = @(t) linspace(leftBound, leftBound+speed*dt, nnodesX);
x_connectivity = [(1:nelsX)', (2:(nelsX+1))'];

% Build global numbering
%numberingSubproblem(i) gives global numbering of node i of Subproblem
numberingX  = 1:nnodesX;
numberingXi = (nnodesX+1):(nnodesX+nnodesXi);

%% Initializing solution
Uxi = icXi(xi_nodesFun(t), 0)';
Ux = icX(x_nodesFun(t), 0)';
U = [Ux; Uxi];
massX = sparse(nnodesX, nnodesX);
massXi = sparse(nnodesXi, nnodesXi);
diffusionX = sparse(nnodesX, nnodesX);
diffusionXi = sparse(nnodesXi, nnodesXi);
advectionXi = sparse(nnodesXi, nnodesXi);
pulse = zeros([nnodesXi, 1]);

xipos = xi_nodesFun(t);
xpos  = x_nodesFun(t);
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
    rloc = h*powerDensity(xiposloc, power, efficiency, radius);%closed integration
    massXi(inodes, inodes) = ...
        massXi(inodes, inodes) + Mloc;
    advectionXi(inodes, inodes) = ...
        advectionXi(inodes, inodes) + Aloc;
    diffusionXi(inodes, inodes) = ...
        diffusionXi(inodes, inodes) + Kloc;
    pulse(inodes) = pulse(inodes) + rloc';
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
        rloc = h*powerDensity(xiposloc, power, efficiency, radius);
        pulse(inodes) = pulse(inodes) + rloc';
    end
    % Assemble subproblems into system
    rhs(numberingX) = rhs(numberingX) + rho*cp*( massX*Ux / dt );
    rhs(numberingXi) = rhs(numberingXi) + rho*cp*( massXi*Uxi / dt ) + pulse;
    lhs(numberingX, numberingX) = lhs(numberingX, numberingX) + rho*cp*(massX/dt) + k*diffusionX;
    lhs(numberingXi, numberingXi) = lhs(numberingXi, numberingXi) + rho*cp*(massXi/dt - speed*advectionXi) + k*diffusionXi;
    % Assemble Dirichlet condition interface
    inodeDirichletXi = xi_connectivity( 1, 1);
    inodeDirichletX = x_connectivity( nelsX, 2 );
    lhs(numberingXi(inodeDirichletXi), :) = 0.0;
    lhs(numberingXi(inodeDirichletXi), numberingX(inodeDirichletX)) = -1.0;
    lhs(numberingXi(inodeDirichletXi), numberingXi(inodeDirichletXi)) = 1.0;
    rhs(numberingXi(inodeDirichletXi)) = 0.0;
    % Assemble Neumann condition interface
    inodeBounX = x_connectivity( nelsX, 2);
    inodesXi = xi_connectivity( 1, :);
    Kboun = [1/2, -1/2];
    lhs(numberingX(inodeBounX), numberingXi(inodesXi)) = ...
        lhs(numberingX(inodeBounX), numberingXi(inodesXi)) + k*Kboun;
    % TODO
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
    plot(posSol+speed*t, sol);
    % xlim([leftBound-speed*(Tfinal+3*dt), rightBound+speed*dt]);
    xlim([leftBound, rightBound]);
    ylim([minY-0.1*range, maxY+0.1*range]);
    pause(0.1);
end

%% HEAT SOURCE
function [pd] = powerDensity(xi, power, efficiency, radius)
    xi0 = 0.0;
    if (abs(xi-xi0)>3*radius)
        pd = 0.0;
    else
        pd = 2*(power*efficiency) / pi / pow2(radius) * exp( - 2*pow2(xi - xi0)/pow2(radius));
    end
end