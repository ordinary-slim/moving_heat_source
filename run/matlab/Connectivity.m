classdef Connectivity < handle
    properties
        odim = -1;
        tdim = -1;
        onents = -1;
        tnents = -1;
        array = [];
    end
    methods
        function [obj] = Connectivity(odim, tdim, onents, tnents, array)
          obj.odim = odim;
          obj.tdim = tdim;
          obj.onents = onents;
          obj.tnents = tnents;
          obj.array = array;
        end
        function [locCon] = getLocalCon(obj, ient)
          locCon = obj.array(ient, find(obj.array(ient, :) >= 0 ) );
        end
        function [transposedCon] = transpose(obj)
            auxcon = cell([obj.tnents, 1]);
            for oent=1:obj.onents
                V = obj.getLocalCon( oent );
                for v=V
                    auxcon{v} = unique([auxcon{v}, oent]);
                end
            end
            transposedCon = Connectivity(obj.tdim, obj.odim, ...
                obj.tnents, obj.onents, aux_array2array( auxcon ) );
        end
    end
end
function [array] = aux_array2array( aux_array )
    maxIncidence = 0;
    for idx=1:length(aux_array)
        if (length(aux_array{idx})>maxIncidence)
            maxIncidence = length(aux_array{idx});
        end
    end
    array = -ones([length(aux_array), maxIncidence]);
    for idx=1:length(aux_array)
        array(idx, 1:length(aux_array{idx})) = aux_array{idx};
    end
end