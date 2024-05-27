classdef HalfHalfScheme < Scheme
    properties
        h
        pulse
        pos1
        pos2
        nnodes1
        nnodes2
        nels1
        nels2
        connectivity1
        connectivity2
        numbering1 
        numbering2
        U1
        U2
        mass1
        diffusion1
        pulse1
        mass2
        diffusion2
        pulse2
        xInterface
        pos%post
        Upos%post
        adimDt
    end
    methods
        function obj = HalfHalfScheme(S)
            obj = obj.load(S);
            obj = obj.initialize(S.icX);
        end
        function obj = initialize(obj, icX)
            obj.h = 1/obj.meshDensity;
            % Initialize meshes
            obj.pos1 = obj.meshLeft();
            obj.pos2  = obj.meshRight();

            obj.nnodes1 = size(obj.pos1,2);
            obj.nnodes2  = size(obj.pos2, 2);
            obj.nels1 = obj.nnodes1-1;
            obj.nels2 = obj.nnodes2-1;
            obj.connectivity1 = [(1:obj.nels1)', (2:(obj.nels1+1))'];
            obj.connectivity2 = [(1:obj.nels2)', (2:(obj.nels2+1))'];

            obj.nnodes = obj.nnodes1 + obj.nnodes2;
            obj.nels   = obj.nels1 + obj.nels2;

            % Build global numbering
            %numberingSubproblem(i) gives global numbering of node i of Subproblem
            obj.numbering1  = 1:obj.nnodes1;
            obj.numbering2  = (obj.nnodes1+1):1:(obj.nnodes1+obj.nnodes2);

            % Initializing matrices
            obj.U1 = icX(obj.pos1, 0)';
            obj.U2 = icX(obj.pos2, 0)';
            obj.U = [obj.U1; obj.U2];
            obj.mass2 = sparse(obj.nnodes2, obj.nnodes2);
            obj.diffusion2 = sparse(obj.nnodes2, obj.nnodes2);
            obj.pulse2 = zeros([obj.nnodes2, 1]);
            obj.mass1 = sparse(obj.nnodes1, obj.nnodes1);
            obj.diffusion1 = sparse(obj.nnodes1, obj.nnodes1);
            obj.pulse1 = zeros([obj.nnodes1, 1]);
        end
        function obj = preLoopAssembly(obj)
          %% Time independent assembly
          % LEFT SUBPROBLEM
          % Mass, diffusion
          for iel=1:obj.nels1
              inodes = obj.connectivity1(iel, :);
              pos1loc = obj.pos1(inodes);
              h = pos1loc(2) - pos1loc(1);
              Mloc = [h/3, h/6; h/6, h/3];
              Kloc = [1/h, -1/h; -1/h, 1/h];
              obj.mass1(inodes, inodes) = ...
                  obj.mass1(inodes, inodes) + Mloc;
              obj.diffusion1(inodes, inodes) = ...
                  obj.diffusion1(inodes, inodes) + Kloc;
          end
          % RIGHT SUBPROBLEM
          % Mass, diffusion
          for iel=1:obj.nels2
              inodes = obj.connectivity2(iel, :);
              pos2loc = obj.pos2(inodes);
              h = pos2loc(2) - pos2loc(1);
              Mloc = [h/3, h/6; h/6, h/3];
              Kloc = [1/h, -1/h; -1/h, 1/h];
              obj.mass2(inodes, inodes) = ...
                  obj.mass2(inodes, inodes) + Mloc;
              obj.diffusion2(inodes, inodes) = ...
                  obj.diffusion2(inodes, inodes) + Kloc;
          end
        end
        function obj = iterate(obj)
          % PRE-ITERATION
          obj.iter = obj.iter + 1;
          "iter " + obj.iter
          obj.t = obj.t + obj.dt;
          obj.x0 = obj.x0 + obj.dt*obj.speed;

          % Reset linear system
          obj.lhs = sparse(obj.nnodes, obj.nnodes);
          obj.rhs = zeros([obj.nnodes, 1]);
          obj.pulse1 = zeros([obj.nnodes1, 1]);
          obj.pulse2 = zeros([obj.nnodes2, 1]);

          %% Time-dependent assembly    
          % Pulse, LEFT
          for iel=1:obj.nels1
              inodes = obj.connectivity1(iel, :);
              pos1loc = obj.pos1(inodes);
              h = pos1loc(2) - pos1loc(1);
              rloc = obj.h*obj.powerDensity(pos1loc, obj.x0);
              obj.pulse1(inodes) = obj.pulse1(inodes) + rloc';
          end
          % Pulse, RIGHT
          for iel=1:obj.nels2
              inodes = obj.connectivity2(iel, :);
              pos2loc = obj.pos2(inodes);
              h = pos2loc(2) - pos2loc(1);
              rloc = h*obj.powerDensity(pos2loc, obj.x0);
              obj.pulse2(inodes) = obj.pulse2(inodes) + rloc';
          end
          % Assemble subproblems into system
          % rhs
          obj.rhs(obj.numbering1) = obj.rhs(obj.numbering1) + obj.rho*obj.cp*( obj.mass1*obj.U1 / obj.dt ) + obj.pulse1;
          obj.rhs(obj.numbering2) = obj.rhs(obj.numbering2) + obj.rho*obj.cp*( obj.mass2*obj.U2 / obj.dt ) + obj.pulse2;
          % lhs
          obj.lhs(obj.numbering1, obj.numbering1) = obj.lhs(obj.numbering1, obj.numbering1) + ...
            obj.rho*obj.cp*(obj.mass1/obj.dt) + obj.k*obj.diffusion1;
          obj.lhs(obj.numbering2, obj.numbering2) = obj.lhs(obj.numbering2, obj.numbering2) + ...
            obj.rho*obj.cp*(obj.mass2/obj.dt) + obj.k*obj.diffusion2;
            
          % EXTERNAL BCs
          % Neumann condition left
          obj.rhs(1) = obj.rhs(1) + obj.k*obj.neumannFluxLeft;
          % Neumann condition right
          obj.rhs(end) = obj.rhs(end) + obj.k*obj.neumannFluxRight; 

          %% INTERFACE
          % Assemble Neumann condition interface RIGHT
          indexNeumannRight = obj.connectivity2( 1, 1);
          indicesNeumannLeft = obj.connectivity1( obj.nels1, :);
          posLeftloc = obj.pos1(indicesNeumannLeft);
          indicesNeumannLeftGlobal = obj.numbering1( indicesNeumannLeft );
          indexNeumannRightGlobal  = obj.numbering2( indexNeumannRight );
          h = pos2loc(2) - pos2loc(1);
          Kboun = -[-1/h, 1/h];
          obj.lhs(indexNeumannRightGlobal, indicesNeumannLeftGlobal) = ...
              obj.lhs(indexNeumannRightGlobal, indicesNeumannLeftGlobal) - obj.k*Kboun;

          % Assemble Dirichlet interface LEFT
          inodeDirichletLeftGlobal = obj.numbering1( obj.connectivity1( obj.nels1, 2) );
          inodeDirichletRightGlobal  = obj.numbering2( obj.connectivity2( 1, 1 ) );
          obj.lhs(inodeDirichletLeftGlobal, :) = 0.0;
          obj.lhs(inodeDirichletLeftGlobal, inodeDirichletLeftGlobal) = 1.0;
          obj.lhs(inodeDirichletLeftGlobal, inodeDirichletRightGlobal) = -1.0;
          obj.rhs(inodeDirichletLeftGlobal) = 0.0;
          %% Solve
          obj.U = obj.lhs \ obj.rhs;
          % POST-ITERATION
          obj.pos = [obj.pos1(1:end-1), obj.pos2];
          obj.Upos = [obj.U(obj.numbering1(1:end-1)); obj.U(obj.numbering2)];
          obj.U1 = obj.U(obj.numbering1);
          obj.U2 = obj.U(obj.numbering2);
          obj.adimDt = obj.speed*obj.dt / obj.radius;
          obj.xInterface = obj.pos1( obj.connectivity1( obj.nels1, 2 ) );%for post
        end
        function [nodes] = meshLeft(obj)
            nodes = (obj.leftBound):obj.h:(obj.leftBound + (obj.rightBound - obj.leftBound)/2);
        end
        function [nodes] = meshRight(obj)
            nodes = (obj.leftBound + (obj.rightBound - obj.leftBound)/2):obj.h:obj.rightBound;
        end
    end
end
