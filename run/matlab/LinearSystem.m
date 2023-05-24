classdef LinearSystem  < handle
    properties
        ndofs
        lhs
        rhs
        U
    end
    methods
        function [obj] = LinearSystem(ndofs)
            obj.ndofs = ndofs;
            obj.cleanup();
        end
        function [obj] = cleanup(obj)
            obj.lhs = sparse(obj.ndofs, obj.ndofs);
            obj.rhs = zeros([obj.ndofs, 1]);
        end
        function [obj] = solve(obj)
            obj.U = obj.lhs \ obj.rhs;
        end
    end
end