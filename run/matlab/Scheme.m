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
      lhs
      rhs
    end
    methods
        function obj = load(obj, workspace)
            load(workspace, "leftBound", "rightBound", "power", ...
                "efficiency", "radius", "cutoffRadius", "x0", ...
                "speed", "rho", "cp", "k", "dt", "meshDensity", "Tfinal");
            % Load problem parameters from a previous MATLAB workspace
            obj.power = power;
            obj.efficiency = efficiency;
            obj.radius = radius;
            obj.x0ini = x0;
            obj.speed = speed;
            obj.rho = rho;
            obj.cp = cp;
            obj.k = k;
            obj.dt = dt;
            obj.Tfinal = Tfinal;
            obj.leftBound = leftBound;
            obj.rightBound = rightBound;
            obj.cutoffRadius = cutoffRadius;
            obj.L = rightBound - leftBound;
            obj.meshDensity = meshDensity;
            obj.x0 = obj.x0ini;
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
