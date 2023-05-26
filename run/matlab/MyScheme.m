classdef MyScheme < Scheme
    properties
        h
        pulse
        xi0
        xipos
        xpos
        nnodesXi
        nnodesX
        nelsXi
        nelsX
        xi_connectivity
        x_connectivity
        numberingX  
        numberingXi 
        Uxi
        Ux
        massX
        diffusionX
        pulseX
        massXi
        diffusionXi
        advectionXi
        pulseXi
        Upos
        pos
        adimDt
        xInterface
        stabilizationLhsXi
        stabilizationRhsXi
        isStabilized = false
        scA = 2.0
        scAD = 4.0
        minLengthSubdomainX = 0.0;
        preAssembled = false;
    end
    methods
        function obj = MyScheme(S)
            obj = obj.load(S);
            if isfield(S, "isStabilized")
                obj.isStabilized = S.isStabilized;
            end
            if isfield(S, "minLengthSubdomainX")
                obj.minLengthSubdomainX = S.minLengthSubdomainX;
            end
            obj = obj.initialize(S.icX, S.icXi);
        end
        function obj = initialize(obj, icX, icXi)
            obj.h = 1/obj.meshDensity;
            obj.xi0 = obj.x0;
            % Initialize meshes
            obj.xipos = obj.getXiNodes();
            obj.xpos  = obj.getXNodes();

            obj.nnodesXi = size(obj.xipos,2);
            obj.nnodesX  = size(obj.xpos, 2);
            obj.nelsXi = obj.nnodesXi-1;
            obj.nelsX = obj.nnodesX-1;
            obj.xi_connectivity = [(1:obj.nelsXi)', (2:(obj.nelsXi+1))'];
            obj.x_connectivity = [(1:obj.nelsX)', (2:(obj.nelsX+1))'];

            obj.nnodes = obj.nnodesXi + obj.nnodesX;
            obj.nels   = obj.nelsXi + obj.nelsX;

            % Build global numbering
            %numberingSubproblem(i) gives global numbering of node i of Subproblem
            obj.numberingX  = 1:obj.nnodesX;
            obj.numberingXi = (obj.nnodesX+1):(obj.nnodesX+obj.nnodesXi);

            % Initializing solution
            obj.Uxi = icXi(obj.xipos, 0)';
            obj.Ux = icX(obj.xpos, 0)';
            obj.U = [obj.Ux; obj.Uxi];
            obj.massX = sparse(obj.nnodesX, obj.nnodesX);
            obj.diffusionX = sparse(obj.nnodesX, obj.nnodesX);
            obj.pulseX = zeros([obj.nnodesX, 1]);
            obj.massXi = sparse(obj.nnodesXi, obj.nnodesXi);
            obj.diffusionXi = sparse(obj.nnodesXi, obj.nnodesXi);
            obj.advectionXi = sparse(obj.nnodesXi, obj.nnodesXi);
            obj.pulseXi = zeros([obj.nnodesXi, 1]);
            % Stabilization
            obj.stabilizationLhsXi = sparse(obj.nnodesXi, obj.nnodesXi);
            obj.stabilizationRhsXi = zeros([obj.nnodesXi, 1]);

        end
        function obj = preLoopAssembly(obj)
          %% Time independent assembly
          % FIXED SUBPROBLEM
          % Assemble mass matrix
          for iel=1:obj.nelsX
              inodes = obj.x_connectivity(iel, :);
              xposloc = obj.xpos(inodes);
              h = xposloc(2) - xposloc(1);
              Mloc = [h/3, h/6; h/6, h/3];
              Kloc = [1/h, -1/h; -1/h, 1/h];
              obj.massX(inodes, inodes) = ...
                  obj.massX(inodes, inodes) + Mloc;
              obj.diffusionX(inodes, inodes) = ...
                  obj.diffusionX(inodes, inodes) + Kloc;
          end
          % XI SUBPROBLEM
          % Assemble mass, advection and diffusion
          for iel=1:obj.nelsXi
              inodes = obj.xi_connectivity(iel, :);
              xiposloc = obj.xipos(inodes);
              h = xiposloc(2) - xiposloc(1);
              Mloc = [h/3, h/6; h/6, h/3];
              Aloc = [-1/2, 1/2; -1/2, 1/2];
              Kloc = [1/h, -1/h; -1/h, 1/h];
              obj.massXi(inodes, inodes) = ...
                  obj.massXi(inodes, inodes) + Mloc;
              obj.advectionXi(inodes, inodes) = ...
                  obj.advectionXi(inodes, inodes) + Aloc;
              obj.diffusionXi(inodes, inodes) = ...
                  obj.diffusionXi(inodes, inodes) + Kloc;
          end
          % Assemble stabilization LHS
          for iel=1:obj.nelsXi
              inodes = obj.xi_connectivity(iel, :);
              xiposloc = obj.xipos(inodes);
              h = xiposloc(2) - xiposloc(1);
              tau = obj.getTau( h );
              Slhsloc = obj.rho * obj.cp * tau * obj.speed^2 / h * [1, -1; -1, 1];
              Srhsloc = tau * (obj.powerDensity(xiposloc, obj.xi0) .* [-1, 1]) * obj.speed;
              obj.stabilizationLhsXi(inodes, inodes) = ...
                  obj.stabilizationLhsXi(inodes, inodes) + Slhsloc;
              obj.stabilizationRhsXi(inodes) = ...
                  obj.stabilizationRhsXi(inodes) + Srhsloc';
          end
        end
        function obj = iterate(obj)
          if ~(obj.preAssembled)
              obj.preLoopAssembly;
              obj.preAssembled = true;
          end
          obj.iter = obj.iter + 1;
          "iter " + obj.iter
          obj.xipos = obj.getXiNodes();
          obj.xpos = obj.getXNodes();
          obj.t = obj.t + obj.dt;
          obj.x0 = obj.x0 + obj.dt*obj.speed;
          obj.lhs = sparse(obj.nnodes, obj.nnodes);
          obj.rhs = zeros([obj.nnodes, 1]);
          obj.pulseXi = zeros([obj.nnodesXi, 1]);
          obj.pulseX = zeros([obj.nnodesX, 1]);
          obj.stabilizationRhsXi = zeros([obj.nnodesXi, 1]);
          %% Time-dependent assembly    
          % Pulse, X
          for iel=1:obj.nelsX
              inodes = obj.x_connectivity(iel, :);
              xposloc = obj.xpos(inodes);
              h = xposloc(2) - xposloc(1);
              rloc = 0.5*h*obj.powerDensity(xposloc, obj.x0);
              obj.pulseX(inodes) = obj.pulseX(inodes) + rloc';
          end
          % Pulse, Xi
          for iel=1:obj.nelsXi
              inodes = obj.xi_connectivity(iel, :);
              xiposloc = obj.xipos(inodes);
              h = xiposloc(2) - xiposloc(1);
              rloc = 0.5*obj.h*obj.powerDensity(xiposloc, obj.xi0);
              obj.pulseXi(inodes) = obj.pulseXi(inodes) + rloc';
          end
          % Stabilization RHS, Xi
          for iel=1:obj.nelsXi
              inodes = obj.xi_connectivity(iel, :);
              xiposloc = obj.xipos(inodes);
              h = xiposloc(2) - xiposloc(1);
              tau = obj.getTau( h );
              Srhsloc = tau * (obj.powerDensity(xiposloc, obj.xi0) .* [-1, 1]) * obj.speed;
              obj.stabilizationRhsXi(inodes) = ...
                  obj.stabilizationRhsXi(inodes) + Srhsloc';
          end
          % Assemble subproblems into system
          obj.rhs(obj.numberingX) = obj.rhs(obj.numberingX) + obj.rho*obj.cp*( obj.massX*obj.Ux / obj.dt ) + obj.pulseX;
          obj.rhs(obj.numberingXi)= obj.rhs(obj.numberingXi) + obj.rho*obj.cp*( obj.massXi*obj.Uxi / obj.dt ) + obj.pulseXi;
          obj.lhs(obj.numberingX, obj.numberingX) = obj.lhs(obj.numberingX, obj.numberingX) + obj.rho*obj.cp*(obj.massX/obj.dt) + obj.k*obj.diffusionX;
          obj.lhs(obj.numberingXi, obj.numberingXi) = obj.lhs(obj.numberingXi, obj.numberingXi) + obj.rho*obj.cp*(obj.massXi/obj.dt - obj.speed*obj.advectionXi) + obj.k*obj.diffusionXi;
            
          if obj.isStabilized
              obj.rhs(obj.numberingXi) = obj.rhs(obj.numberingXi) + obj.stabilizationRhsXi;
              obj.lhs(obj.numberingXi, obj.numberingXi) = obj.lhs(obj.numberingXi, obj.numberingXi) + obj.stabilizationLhsXi;
          end
          

          % EXTERNAL BCs
          % Neumann condition left
          obj.rhs(1) = obj.rhs(1) + obj.k*obj.neumannFluxLeft;
          % Neumann condition right
          obj.rhs(end) = obj.rhs(end) + obj.k*obj.neumannFluxRight; 

          %% INTERFACE
          % Assemble Neumann condition interface LEFT
          % inodeBounX = obj.x_connectivity( obj.nelsX, 2);
          % inodesXi = obj.xi_connectivity( 1, :);
          % xiposloc = obj.xipos(inodesXi);
          % h = xiposloc(2) - xiposloc(1);
          % Kboun = [-1/h, 1/h];
          % obj.lhs(obj.numberingX(inodeBounX), obj.numberingXi(inodesXi)) = ...
          %     obj.lhs(obj.numberingX(inodeBounX), obj.numberingXi(inodesXi)) - obj.k*Kboun;

          % Assemble Dirichlet interface RIGHT
          % inodeDirichletXi = obj.xi_connectivity( 1, 1);
          % inodeDirichletX = obj.x_connectivity( obj.nelsX, 2 );
          % obj.lhs(obj.numberingXi(inodeDirichletXi), :) = 0.0;
          % obj.lhs(obj.numberingXi(inodeDirichletXi), obj.numberingX(inodeDirichletX)) = -1.0;
          % obj.lhs(obj.numberingXi(inodeDirichletXi), obj.numberingXi(inodeDirichletXi)) = 1.0;
          % obj.rhs(obj.numberingXi(inodeDirichletXi)) = 0.0;

          % Assemble Neumann condition interface RIGHT
          inodeBounXi_Neumann = obj.xi_connectivity( 1, 1);
          inodesX_Neumann = obj.x_connectivity( obj.nelsX, :);
          xposloc = obj.xpos(inodesX_Neumann);
          inodeBounXi_Neumann_global = obj.numberingXi( inodeBounXi_Neumann );
          inodesX_Neumann_global = obj.numberingX( inodesX_Neumann );
          h = xposloc(2) - xposloc(1);
          Kboun = -[-1/h, 1/h];
          obj.lhs(inodeBounXi_Neumann_global, inodesX_Neumann_global) = ...
              obj.lhs(inodeBounXi_Neumann_global, inodesX_Neumann_global) - obj.k*Kboun;
          % TODO, fix here
          % Assemble Dirichlet interface LEFT
          inodeDirichletXi_global = obj.numberingXi( obj.xi_connectivity( 1, 1) );
          inodeDirichletX_global = obj.numberingX( obj.x_connectivity( obj.nelsX, 2 ) );
          obj.lhs(inodeDirichletX_global, :) = 0.0;
          obj.lhs(inodeDirichletX_global, inodeDirichletX_global) = 1.0;
          obj.lhs(inodeDirichletX_global, inodeDirichletXi_global) = -1.0;
          obj.rhs(inodeDirichletX_global) = 0.0;
          %% Solve
          obj.U = obj.lhs \ obj.rhs;
          % POST-ITERATION
          obj.pos = [obj.xpos(1:end-1)-obj.speed*obj.t, obj.xipos];
          obj.Upos = [obj.U(obj.numberingX(1:end-1)); obj.U(obj.numberingXi)];
          obj.Ux = obj.U(obj.numberingX);
          obj.Uxi = interp1(obj.pos, obj.Upos, obj.getXiNodes(), "linear", "extrap")';%gotta work on this
          obj.adimDt = obj.speed*obj.dt / obj.radius;
          obj.xInterface = obj.xpos( obj.x_connectivity( obj.nelsX, 2 ) );%for post
        end
        function [tau] = getTau(obj, h)
          advectionEstimate = h / obj.scA / (obj.rho*obj.cp*abs(obj.speed));
          if obj.k>0
            diffusionEstimate = h.^2 / obj.scAD / obj.k;
            tau = 1 / (1/advectionEstimate + 1/diffusionEstimate);
          else
            tau = advectionEstimate;
          end
        end
        function lengthSubdomainX = getLengthSubdomainX(obj)
            lengthSubdomainX = (obj.speed*obj.dt);
            lengthSubdomainX = max( obj.minLengthSubdomainX, lengthSubdomainX);
        end
        function [nodes] = getXiNodes(obj)   
            lengthSubdomainX = obj.getLengthSubdomainX;
            leftBoundXi = round((-obj.speed*(obj.t+obj.dt)+obj.leftBound+lengthSubdomainX),5);
            nodes = leftBoundXi:obj.h:(-obj.speed*(obj.t+obj.dt) + obj.rightBound);
        end
        function [nodes] = getXNodes(obj)
            lengthSubdomainX = obj.getLengthSubdomainX;
            nodes = (obj.leftBound):obj.h:(obj.leftBound+lengthSubdomainX);
        end
    end
end
