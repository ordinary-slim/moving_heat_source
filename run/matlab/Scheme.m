classdef Scheme < handle
    properties
      power = 0.0;
      efficiency = 1.0;
      radius = 1.0;
      x0ini = 0.0;
      x0
      speed = 0.0;
      rho = 1.0;
      cp = 1.0;
      k = 1.0;
      dt = 0.1;
      Tfinal = 1.0;
      leftBound = 0.0;
      rightBound = 1.0;
      L = 1.0;
      cutoffRadius = 3;
      meshDensity = 1;
      t = 0.0;
      nnodes = 2;
      nels = 1;
      U = [];
      iter = 0;
      neumannFluxLeft = 0.0;
      neumannFluxRight = 0.0;
      lhs
      rhs
    end
    methods
        function obj = load(obj, S)
            % Load problem parameters from a previous MATLAB workspace
            obj.power = S.power;
            obj.efficiency = S.efficiency;
            obj.radius = S.radius;
            obj.x0ini = S.x0;
            obj.speed = S.speed;
            obj.rho = S.rho;
            obj.cp = S.cp;
            obj.k = S.k;
            obj.dt = S.dt;
            obj.Tfinal = S.Tfinal;
            obj.leftBound = S.leftBound;
            obj.rightBound = S.rightBound;
            obj.cutoffRadius = S.cutoffRadius;
            obj.L = S.rightBound - S.leftBound;
            obj.meshDensity = S.meshDensity;
            obj.x0 = obj.x0ini;
            if isfield(S, "neumannFluxLeft")
                obj.neumannFluxLeft = S.neumannFluxLeft;
            end
            if isfield(S, "neumannFluxRight")
                obj.neumannFluxRight = S.neumannFluxRight;
            end
        end
        function [pd] = powerDensity(obj, x, x0)
            if (abs(x-x0)>obj.cutoffRadius)
                pd = 0.0;
            else
                pd = 2*(obj.power*obj.efficiency) / pi / obj.radius^2 * exp( - 2*(x - x0).^2/obj.radius^2);
            end
        end
    end
end
