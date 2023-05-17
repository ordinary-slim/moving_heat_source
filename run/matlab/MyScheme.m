classdef MyScheme < Scheme
    properties
        h
        connectivity
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
    end
    methods
        function obj = MyScheme(workspace)
            obj = obj.load(workspace);
            load(workspace, "icX", "icXi");
            obj = obj.initialize(icX, icXi);
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
          % Assemble advection and mass mat
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
        end
        function obj = iterate(obj)
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
          %% Time-dependent assembly    
          % Pulse, X
          for iel=1:obj.nelsX
              inodes = obj.x_connectivity(iel, :);
              xposloc = obj.xpos(inodes);
              h = xposloc(2) - xposloc(1);
              rloc = h*obj.powerDensity(xposloc, obj.x0);
              obj.pulseX(inodes) = obj.pulseX(inodes) + rloc';
          end
          % Pulse, Xi
          for iel=1:obj.nelsXi
              inodes = obj.xi_connectivity(iel, :);
              xiposloc = obj.xipos(inodes);
              h = xiposloc(2) - xiposloc(1);
              rloc = obj.h*obj.powerDensity(xiposloc, obj.xi0);
              obj.pulseXi(inodes) = obj.pulseXi(inodes) + rloc';
          end
          % Assemble subproblems into system
          obj.rhs(obj.numberingX) = obj.rhs(obj.numberingX) + obj.rho*obj.cp*( obj.massX*obj.Ux / obj.dt ) + obj.pulseX;
          obj.rhs(obj.numberingXi)= obj.rhs(obj.numberingXi) + obj.rho*obj.cp*( obj.massXi*obj.Uxi / obj.dt ) + obj.pulseXi;
          obj.lhs(obj.numberingX, obj.numberingX) = obj.lhs(obj.numberingX, obj.numberingX) + obj.rho*obj.cp*(obj.massX/obj.dt) + obj.k*obj.diffusionX;
          obj.lhs(obj.numberingXi, obj.numberingXi) = obj.lhs(obj.numberingXi, obj.numberingXi) + obj.rho*obj.cp*(obj.massXi/obj.dt - obj.speed*obj.advectionXi) + obj.k*obj.diffusionXi;


          %% INTERFACE
          % Assemble Neumann condition interface LEFT
          % inodeBounX = x_connectivity( nelsX, 2);
          % inodesXi = xi_connectivity( 1, :);
          % xiposloc = xipos(inodesXi);
          % h = xiposloc(2) - xiposloc(1);
          % Kboun = [-1/h, 1/h];
          % lhs(numberingX(inodeBounX), numberingXi(inodesXi)) = ...
          %     lhs(numberingX(inodeBounX), numberingXi(inodesXi)) - k*Kboun;

          % Assemble Dirichlet interface RIGHT
          % inodeDirichletXi = xi_connectivity( 1, 1);
          % inodeDirichletX = x_connectivity( nelsX, 2 );
          % lhs(numberingXi(inodeDirichletXi), :) = 0.0;
          % lhs(numberingXi(inodeDirichletXi), numberingX(inodeDirichletX)) = -1.0;
          % lhs(numberingXi(inodeDirichletXi), numberingXi(inodeDirichletXi)) = 1.0;
          % rhs(numberingXi(inodeDirichletXi)) = 0.0;

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
          obj.Uxi = interp1(obj.pos, obj.Upos, obj.getXiNodes())';%gotta work on this
          obj.adimDt = obj.speed*obj.dt / obj.radius;
        end
        function [nodes] = getXiNodes(obj)
            nodes = (-obj.speed*obj.t+obj.leftBound):obj.h:(-obj.speed*(obj.t+obj.dt) + obj.rightBound);
        end
        function [nodes] = getXNodes(obj)
            nodes = (obj.leftBound):obj.h:(obj.leftBound+obj.speed*obj.dt);
        end
    end
end
