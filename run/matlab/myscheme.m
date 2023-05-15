clc; clear all; close all;
%% PARAMS
L = 1;
V = 0.5;
leftBound = 0;
rightBound = leftBound + L;
% IC
m = 11;%slope of initial condition
c = 0;
exactSolXi = @(xi, t) c + m*(xi + V*t);
exactSolX = @(x, t) c + m*(x);
% exactSolXi = @(xi, t) max(0, 1 - abs(xi+V*t));
% exactSolX = @(x, t) max(0, 1 - abs(xi+V*t));
% MATERIAL
rho = 1;
cp = 1;
k = 0;

%% DISCRETIZATION
nelsXi = 8;%# elements
nnodesXi = nelsXi + 1;
nelsX = 2;%# elements
nnodesX = nelsX + 1;
nnodes = nnodesXi + nnodesX;
nels   = nelsXi + nelsX;
dt = 0.1;
dx = L/nels;
tol = 1e-7;
t = 0.0;
Tfinal = 1;

%% MESHING
xi_nodesFun = @(t) linspace(-V*t + leftBound,-V*(t+dt) + rightBound, nnodesXi);
xi_connectivity = [(1:nelsXi)', (2:(nelsXi+1))'];
x_nodesFun = @(t) linspace(leftBound, V*dt, nnodesX);
x_connectivity = [(1:nelsX)', (2:(nelsX+1))'];

% Build global numbering
%numberingSubproblem(i) gives global numbering of node i of Subproblem
numberingX  = 1:nnodesX;
numberingXi = (nnodesX+1):(nnodesX+nnodesXi);

%% Initializing solution
Uxi = exactSolXi(xi_nodesFun(t), 0)';
Ux = exactSolX(x_nodesFun(t), 0)';
U = [Ux; Uxi];
massX = sparse(nnodesX, nnodesX);
massXi = sparse(nnodesXi, nnodesXi);
diffusionX = sparse(nnodesX, nnodesX);
diffusionXi = sparse(nnodesXi, nnodesXi);
advectionXi = sparse(nnodesXi, nnodesXi);

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
    massXi(inodes, inodes) = ...
        massXi(inodes, inodes) + Mloc;
    advectionXi(inodes, inodes) = ...
        advectionXi(inodes, inodes) + Aloc;
    diffusionXi(inodes, inodes) = ...
        diffusionXi(inodes, inodes) + Kloc;
end

%% Time loop
iter = 0;
while t < Tfinal
    iter = iter + 1;
    "iter " + iter
    xipos = xi_nodesFun(t);
    xpos = x_nodesFun(t);
    t = t + dt;
    %% Assembly    
    lhs = sparse(nnodes, nnodes);
    rhs = zeros([nnodes, 1]);
    % Assemble subproblems into system
    rhs(numberingX) = rhs(numberingX) + rho*cp*( massX*Ux / dt );
    rhs(numberingXi) = rhs(numberingXi) + rho*cp*( massXi*Uxi / dt );
    lhs(numberingX, numberingX) = lhs(numberingX, numberingX) + rho*cp*(massX/dt) + k*diffusionX;
    lhs(numberingXi, numberingXi) = lhs(numberingXi, numberingXi) + rho*cp*(massXi/dt - V*advectionXi) + k*diffusionXi;
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
    lhs(inodeBounX, inodesXi) = lhs(inodeBounX, inodesXi) + k*Kboun;
    % TODO
    %% Solve
    U = lhs \ rhs;

    %% pos-iteration
    posSol = [xpos(1:end-1)-V*t, xipos];
    sol = [U(numberingX(1:end-1)); U(numberingXi)];
    Ux = U(numberingX);
    Uxi = interp1(posSol, sol, xi_nodesFun(t))';%gotta work on this
    Uxi = exactSolXi(xi_nodesFun(t), t)';

    
    plot(posSol, sol);
    xlim([-0.75*L, 1.25*L]);
    ylim([-0.1*m, 1.1*L*m]);
    pause(0.1);
end