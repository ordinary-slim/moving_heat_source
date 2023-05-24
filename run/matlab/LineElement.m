classdef LineElement  < handle
    properties
        pos = zeros([2, 1]);
        gpos = zeros([2, 1]);
        gpweight = [0.5, 0.5];%closed integration
        con = zeros([2, 1]);
        vol = 0;
        nnodes = 2;
        ngpoints = 2;
        baseFuns = {};
        gradBaseFuns = {};
        centroid = 0.0;
    end
    methods
        function [obj] = LineElement( pos, con )
            obj.pos = pos;
            obj.gpos = pos;%closed integration
            obj.con = con;
            obj.initialize();
        end
        function [obj] = initialize(obj)
            obj.computeVolume();
            obj.computeCentroid();
            obj.baseFuns = {@(x) 1-(x-obj.pos(1))/obj.vol, @(x) (x-obj.pos(1))/obj.vol};
            obj.gradBaseFuns = {@(x) -1/obj.vol, @(x) 1/obj.vol};
        end
        function [obj] = computeVolume(obj)
            obj.vol = obj.pos(2) - obj.pos(1);
        end
        function [obj] = computeCentroid(obj)
            obj.centroid = (obj.pos(2) + obj.pos(1))/2;
        end
    end
end