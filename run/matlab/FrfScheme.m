classdef FrfScheme < Scheme
    properties
        xpos
        h
        connectivity
        massX = [];
        diffusionX = [];
        pulse
    end
    methods
        function obj = FrfScheme(workspace)
            obj = obj.load(workspace);
            load(workspace, "icX");
            obj = obj.initialize(icX);
        end
        function obj = initialize(obj, icX)
            obj.h = 1/obj.meshDensity;
            obj.xpos  = obj.getNodes();
            obj.nnodes  = size(obj.xpos, 2);
            obj.nels = obj.nnodes-1;
            obj.connectivity = [(1:obj.nels)', (2:(obj.nels+1))'];
            obj.U = icX(obj.xpos, 0)';
            obj.massX = sparse(obj.nnodes, obj.nnodes);
            obj.diffusionX = sparse(obj.nnodes, obj.nnodes);
        end
        function obj = preLoopAssembly(obj)
            for iel=1:obj.nels
                inodes = obj.connectivity(iel, :);
                xposloc = obj.xpos(inodes);
                h = xposloc(2) - xposloc(1);
                Mloc = [h/3, h/6; h/6, h/3];
                Kloc = [1/h, -1/h; -1/h, 1/h];
                obj.massX(inodes, inodes) = ...
                    obj.massX(inodes, inodes) + Mloc;
                obj.diffusionX(inodes, inodes) = ...
                    obj.diffusionX(inodes, inodes) + Kloc;
            end
        end
        function obj = iterate(obj)
                obj.iter = obj.iter + 1;
                "iter " + obj.iter
                obj.t = obj.t + obj.dt;
                obj.x0 = obj.x0 + obj.dt*obj.speed;
                %% Assembly    
                obj.lhs = sparse(obj.nnodes, obj.nnodes);
                obj.rhs = zeros([obj.nnodes, 1]);
                obj.pulse = zeros([obj.nnodes, 1]);
                % Assemble heat source
                for iel=1:obj.nels
                    inodes = obj.connectivity(iel, :);
                    xloc = obj.xpos( inodes );
                    rloc = obj.h*obj.powerDensity(xloc, obj.x0);
                    obj.pulse(inodes) = obj.pulse(inodes) + rloc';
                end
                % Assemble subproblems into system
                obj.rhs = obj.rhs + obj.rho*obj.cp*( obj.massX*obj.U / obj.dt ) + obj.pulse;
                obj.lhs = obj.lhs + obj.rho*obj.cp*(obj.massX/obj.dt) + obj.k*obj.diffusionX;
                %% Solve
                obj.U = obj.lhs \ obj.rhs;
        end
        function [nodes] = getNodes(obj)
            nodes = (obj.leftBound):obj.h:(obj.rightBound);
        end
    end
end
