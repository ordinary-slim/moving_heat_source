classdef FRFScheme < handle
    properties
        ls LinearSystem
        numberingDofs
        mesh Mesh
        problem Subproblem
    end
    methods
        function [obj] = FRFScheme(params)
            %MESHING
            nelsPart = (params.rightBound - params.leftBound) * params.meshDensity;
            nelsHS = params.radiusSubdomain*2*params.meshDensity;
            obj.mesh = Mesh( params.leftBound, params.rightBound, nelsPart );            
            
            ndofs = obj.mesh.nnodes;
            %NUMBERING DOFS
            obj.numberingDofs = 1:obj.mesh.nnodes;
            obj.ls = LinearSystem( ndofs );
            
            % Subproblems
            params.vadv = 0.0;
            params.vsource = params.speed;
            params.vdomain = 0.0;
            obj.problem = Subproblem(obj.mesh, obj.numberingDofs, params);
        end
        function [obj] = iterate(obj)
          %% pre-iteration
          obj.ls.cleanup();
          obj.problem.preIterate();
          %% Assembly
          obj.problem.assemblePDE(obj.ls);
          %% Solve
          obj.ls.solve();
          %% post-iteration
          obj.problem.U = obj.ls.U(obj.problem.numbering);
        end
        function [t] = getTime(obj)
            t = obj.problem.time;
        end
    end
end
