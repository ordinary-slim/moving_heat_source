classdef Subproblem  < handle
    % Member of a 1D coupled problem
  properties
    mesh Mesh
    numbering = [];
    U = [];
    dt = 0.0;
    time  = 0.0;
    vdomain = 0.0;
    % Material properties
    rho = 0.0;
    cp = 0.0;
    k = 0.0;
    vadv = 0.0;
    % Heat source
    x0 = 0.0;
    power = 0.0;
    radius = 0.0;
    vsource = 0.0;
  end
  methods
      function [obj] = Subproblem( mesh, numbering, params )
        obj.mesh = mesh;
        obj.numbering = numbering;
        obj.U = params.ic( mesh.pos );
        obj.dt = params.dt;
        obj.vdomain = params.vdomain;
        % Material properties
        obj.rho = params.rho;
        obj.cp = params.cp;
        obj.k = params.k;
        obj.vadv = params.vadv;
        % Heat source
        obj.x0 = params.x0;
        obj.power = params.power;
        obj.radius = params.radius;
        obj.vsource = params.vsource;
      end
      function [obj]= preIterate(obj)
          obj.time = obj.time + obj.dt;
          obj.x0 = obj.x0 + obj.dt * obj.vsource;
          obj.mesh.updatePosFixed( obj.mesh.shiftFixed + obj.dt * obj.vdomain);
      end
      function assemblePDE( obj, linearSystem )
          arguments
              obj Subproblem
              linearSystem LinearSystem
          end
          obj.assembleInteriorSpatial(linearSystem);
          obj.assembleTimeDerivative(linearSystem);
      end
      function assembleInteriorSpatial(obj, linearSystem)
          arguments
              obj Subproblem
              linearSystem LinearSystem
          end
          for ielem=1:obj.mesh.nels
              if ~(obj.mesh.activeElements(ielem))
                  continue;
              end
              % Load element
              e = obj.mesh.getElement( ielem );
              % Compute contributions
              Aloc = zeros([e.nnodes, e.nnodes]);
              rloc = zeros([e.nnodes, 1]);
              for igp=1:e.ngpoints
                  xgp = e.gpos(igp);
                  for inode=1:e.nnodes
                      rloc(inode) = rloc(inode) + ...
                          e.baseFuns{inode}(xgp)*obj.powerDensity(xgp)*...
                          e.gpweight(igp)*e.vol;
                      for jnode=1:e.nnodes
                          %DIFFUSION
                          Aloc(inode, jnode) = Aloc(inode, jnode) + ...
                              obj.k*...
                              e.gradBaseFuns{inode}(xgp) * e.gradBaseFuns{jnode}(xgp)*...
                              e.gpweight(igp)*e.vol;
                          %ADVECTION
                          Aloc(inode, jnode) = Aloc(inode, jnode) + ...
                              obj.rho*obj.cp*obj.vadv*...
                              e.baseFuns{inode}(xgp) * e.gradBaseFuns{jnode}(xgp)*...
                              e.gpweight(igp)*e.vol;
                      end
                  end
              end
              % Assemble
              globalCon = obj.numbering(e.con);
              linearSystem.lhs(globalCon, globalCon) = linearSystem.lhs(globalCon, globalCon) + ...
                  Aloc;
              linearSystem.rhs(globalCon) = linearSystem.rhs(globalCon) + ...
                  rloc;
          end
      end
      function assembleTimeDerivative(obj, linearSystem)
          arguments
              obj Subproblem
              linearSystem LinearSystem
          end
          % Time derivative
          linearSystem.lhs(obj.numbering, obj.numbering) = linearSystem.lhs(obj.numbering, obj.numbering) + ...
                  obj.rho*obj.cp*obj.mesh.mass / obj.dt;
          linearSystem.rhs(obj.numbering) = linearSystem.rhs(obj.numbering) + ...
              obj.rho*obj.cp*(obj.mesh.mass * obj.U) / obj.dt;
      end
      function assembleInactiveNodes(obj, linearSystem)
          arguments
              obj Subproblem
              linearSystem LinearSystem
          end
          for inode=1:obj.mesh.nnodes
              if ~(obj.mesh.activeNodes(inode))
                  globalIdx = obj.numbering(inode);
                  linearSystem.lhs(globalIdx, :) = 0;
                  linearSystem.lhs(globalIdx, globalIdx) = 1;
                  linearSystem.rhs(globalIdx) = obj.U(inode);
              end
          end
      end
      function u = evaluate( obj, x )
          ielem = obj.mesh.findOwnerElement( x );
          e = obj.mesh.getElement( ielem );
          coeffs = [e.baseFuns{1}(x), e.baseFuns{2}(x)];
          u = dot( coeffs, obj.U(e.con));
      end
      function [pd] = powerDensity(obj, x)
          pd = 2*(obj.power) / pi / obj.radius^2 * exp( - 2*(x - obj.x0).^2/obj.radius^2);
      end
  end
end
