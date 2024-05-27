clear;
inputdset = "myDDPrototype.mat";
params = load(inputdset, "leftBound", "rightBound", "power", ...
    "efficiency", "radius", "cutoffRadius", "x0", ...
    "speed", "rho", "cp", "k", "dt", "meshDensity", "Tfinal", "icX", "icXi");

% PARAMS
params.Tenv = 25;
params.power = 100.0;
params.meshDensity =2;
params.Tfinal = 3;
% params.k = 0.0;
params.leftBound = -25;
params.rightBound = 25;
params.x0 = -15;
params.radiusSubdomain = 5*params.radius;
setDt( params, 2 );
params.ic = @(x) params.Tenv*ones(size(x));


%MESHING
nelsPart = (params.rightBound - params.leftBound) * params.meshDensity;
meshPart = Mesh( params.leftBound, params.rightBound, nelsPart );
nelsHS = radiusSubdomain*2*params.meshDensity;
meshSubdomain = Mesh( params.x0-params.radiusSubdomain, params.x0+params.radiusSubdomain, nelsHS );


ndofs = meshPart.nnodes + meshSubdomain.nnodes;
%NUMBERING DOFS
numberingPart = 1:meshPart.nnodes;
numberingSubdomain = (numberingPart(end)+1):(numberingPart(end)+meshSubdomain.nnodes);
ls = LinearSystem( ndofs );

% Subproblems
params.ic = @(x) params.Tenv*ones(size(x));
params.vadv = 0.0;
params.vsource = params.speed;
params.vdomain = 0.0;
problemPart = Subproblem(meshPart, numberingPart, params);
params.vadv = -params.speed;
params.vsource = 0.0;
params.vdomain = params.speed;
problemHS = Subproblem(meshSubdomain, numberingSubdomain, params);

tol = 1e-7;
figure
while params.Tfinal-tol > problemPart.time
    %% pre-iteration
    ls.cleanup();
    problemPart.preIterate();
    problemHS.preIterate();
    % Domain OPS
    [~, partGammaNodes] = problemPart.mesh.substract( problemHS.mesh );
    problemHS.mesh.intersect( problemPart.mesh );
    hsGammaNodes = problemHS.mesh.bounNodes;%TODO: fix this
    %% Assembly
    problemPart.assemblePDE(ls);
    problemHS.assemblePDE(ls);
    % INTERFACE
    % problemPart receives DIRICHLET from HS
    for inode=partGammaNodes
        posNode_other = problemPart.mesh.posFixed(inode) - problemHS.mesh.shiftFixed;
        iowner = problemHS.mesh.findOwnerElement( posNode_other );
        e = problemHS.mesh.getElement( iowner );
        coeffs = [e.baseFuns{1}(posNode_other), e.baseFuns{2}(posNode_other)];

        globIdx_receiveDirichlet = problemPart.numbering( inode );
        globIdx_sendDirichlet = problemHS.numbering( e.con );

        ls.rhs(globIdx_receiveDirichlet) = 0;
        ls.lhs(globIdx_receiveDirichlet, :) = 0.0;
        ls.lhs(globIdx_receiveDirichlet, globIdx_receiveDirichlet) = -1;
        ls.lhs(globIdx_receiveDirichlet, globIdx_sendDirichlet) = coeffs;
    end
    % problemHS receives NEUMANN from Part
    for inode=hsGammaNodes
        posNode = problemHS.mesh.pos( inode );
        posNode_other = problemHS.mesh.posFixed(inode) - problemPart.mesh.shiftFixed;
        % Compute normal
        parentEl = meshSubdomain.getElement( meshSubdomain.conPointCell.getLocalCon( inode ) );
        normal = posNode - parentEl.centroid;
        normal = normal / abs(normal);
        % Get corresponding element neighbour mesh
        iowner = problemPart.mesh.findOwnerElement( posNode_other );
        e = problemPart.mesh.getElement( iowner );
        %TODO: THink about conductivities here!
        dn_coeffs = normal*[e.gradBaseFuns{1}(posNode_other), e.gradBaseFuns{2}(posNode_other)];

        globIdx_receiveNeumann = problemHS.numbering( inode );
        globIdx_sendNeumann = problemPart.numbering( e.con );
        
        ls.lhs(globIdx_receiveNeumann, globIdx_sendNeumann) = ls.lhs(globIdx_receiveNeumann, globIdx_sendNeumann) - ...
            problemHS.k * dn_coeffs;
    end

    problemPart.assembleInactiveNodes(ls);
    problemHS.assembleInactiveNodes(ls);
    %% Solve
    ls.solve();
    %% post-iteration
    problemPart.U = ls.U(problemPart.numbering);
    problemHS.U = ls.U(problemHS.numbering);
    % Interpolate inactive nodes of problemPart
    inactiveNodesPart = find(~problemPart.mesh.activeNodes);
    for inode=inactiveNodesPart'
        posNode_other = problemPart.mesh.posFixed(inode) - problemHS.mesh.shiftFixed;
        problemPart.U(inode) =  problemHS.evaluate( posNode_other );
    end
    %% plot
    hold off
    plot(problemPart.mesh.posFixed, problemPart.U, 'r')
    hold on
    pause(0.25)
end


function [dt] = setDt( S, adimR)
    dt = S.radius / S.speed * adimR; 
end