clc; clear all;
%% PARAMS
L = 50;
leftBound = -L/2;
rightBound = leftBound + L;

% HEAT SOURCE
power = 100.0;
efficiency = 1.0;
radius = 1;
x0 = -5.0;
speed = 10;
% IC
icX = @(x, t) 25*ones(size(x));

% MATERIAL
rho = 1;
cp = 1;
k = 1;

%% DISCRETIZATION
% TIME
tol = 1e-7;
t = 0.0;
dt = 0.1;
Tfinal = 2.0;
% SPACE
meshDen = 2;
h = 1/meshDen;
x_nodesFun = @(t) (leftBound):h:(rightBound);
xpos  = x_nodesFun(t);
nnodes  = size(xpos, 2);
nels = nnodes-1;
connectivity = [(1:nels)', (2:(nels+1))'];


%% Initializing solution
U = icX(xpos, 0)';
massX = sparse(nnodes, nnodes);
diffusionX = sparse(nnodes, nnodes);

% FIXED SUBPROBLEM
% Assemble mass matrix
for iel=1:nels
    inodes = connectivity(iel, :);
    xposloc = xpos(inodes);
    h = xposloc(2) - xposloc(1);
    Mloc = [h/3, h/6; h/6, h/3];
    Kloc = [1/h, -1/h; -1/h, 1/h];
    massX(inodes, inodes) = ...
        massX(inodes, inodes) + Mloc;
    diffusionX(inodes, inodes) = ...
        diffusionX(inodes, inodes) + Kloc;
end

%% Time loop
iter = 0;
figure
while t < Tfinal
    iter = iter + 1;
    "iter " + iter
    t = t + dt;
    x0 = x0 + dt*speed;
    %% Assembly    
    lhs = sparse(nnodes, nnodes);
    rhs = zeros([nnodes, 1]);
    pulse = zeros([nnodes, 1]);
    % Assemble heat source
    for iel=1:nels
        inodes = connectivity(iel, :);
        xloc = xpos( inodes );
        rloc = h*powerDensity(xloc, x0, power, efficiency, radius);
        pulse(inodes) = pulse(inodes) + rloc';
    end
    % Assemble subproblems into system
    rhs = rhs + rho*cp*( massX*U / dt ) + pulse;
    lhs = lhs + rho*cp*(massX/dt) + k*diffusionX;
    %% Solve
    U = lhs \ rhs;

    minY = min(U);
    maxY = max(U);
    range = maxY - minY;
    plot(xpos, U);
    % xlim([leftBound-speed*(Tfinal+3*dt), rightBound+speed*dt]);
    xlim([leftBound, rightBound]);
    ylim([minY-0.1*range, maxY+0.1*range]);
    pause(0.05);
end

%% HEAT SOURCE
function [pd] = powerDensity(x, x0, power, efficiency, radius)
    if (abs(x-x0)>3*radius)
        pd = 0.0;
    else
        pd = 2*(power*efficiency) / pi / radius^2 * exp( - 2*(x - x0).^2/radius^2);
    end
end