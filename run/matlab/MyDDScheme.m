classdef MyDDScheme < handle
    properties
        ls LinearSystem
        numberingPart
        numberingSubdomain
        meshPart Mesh
        meshSubdomain Mesh
        problemPart Subproblem
        problemHS Subproblem
        posInterface = [];
    end
    methods
        function [obj] = MyDDScheme(params)
            %MESHING
            nelsPart = (params.rightBound - params.leftBound) * params.meshDensity;
            nelsHS = params.radiusSubdomain*2*params.meshDensity;
            obj.meshPart = Mesh( params.leftBound, params.rightBound, nelsPart );
            obj.meshSubdomain = Mesh( params.x0-params.radiusSubdomain, params.x0+params.radiusSubdomain, nelsHS );
            
            
            ndofs = obj.meshPart.nnodes + obj.meshSubdomain.nnodes;
            %NUMBERING DOFS
            obj.numberingPart = 1:obj.meshPart.nnodes;
            obj.numberingSubdomain = (obj.numberingPart(end)+1):(obj.numberingPart(end)+obj.meshSubdomain.nnodes);
            obj.ls = LinearSystem( ndofs );
            
            % Subproblems
            params.vadv = 0.0;
            params.vsource = params.speed;
            params.vdomain = 0.0;
            obj.problemPart = Subproblem(obj.meshPart, obj.numberingPart, params);
            params.vadv = -params.speed;
            params.vsource = 0.0;
            params.vdomain = params.speed;
            obj.problemHS = Subproblem(obj.meshSubdomain, obj.numberingSubdomain, params);
        end
        function [obj] = iterate(obj)
          %% pre-iteration
          obj.ls.cleanup();
          % Motion is done here
          obj.problemHS.mesh.intersect( obj.problemPart.mesh, true );
          obj.problemPart.preIterate();
          obj.problemHS.preIterate();
          % Domain OPS
          obj.problemHS.mesh.intersect( obj.problemPart.mesh, false );
          obj.problemPart.mesh.substract( obj.problemHS.mesh );
          % Determine Gamma boundary
          partGammaNodes = [];
          for inode=obj.problemPart.mesh.bounNodes
              posNode_other = obj.problemPart.mesh.posFixed(inode) - obj.problemHS.mesh.shiftFixed;
              if obj.meshSubdomain.findOwnerElement(posNode_other) >= 0
                  partGammaNodes = [partGammaNodes, inode];
              end
          end
          hsGammaNodes = [];
          for inode=obj.problemHS.mesh.bounNodes
              if obj.meshPart.findOwnerElement( obj.problemHS.mesh.posFixed(inode) ) >= 0
                  hsGammaNodes = [hsGammaNodes, inode];
              end
          end
          %% Assembly
          obj.problemPart.assemblePDE(obj.ls);
          obj.problemHS.assemblePDE(obj.ls);
          % INTERFACE
          % obj.problemPart receives DIRICHLET from HS
          for inode=partGammaNodes
              posNode_other = obj.problemPart.mesh.posFixed(inode) - obj.problemHS.mesh.shiftFixed;
              iowner = obj.problemHS.mesh.findOwnerElement( posNode_other );
              e = obj.problemHS.mesh.getElement( iowner );
              coeffs = [e.baseFuns{1}(posNode_other), e.baseFuns{2}(posNode_other)];

              globIdx_receiveDirichlet = obj.problemPart.numbering( inode );
              globIdx_sendDirichlet = obj.problemHS.numbering( e.con );

              obj.ls.rhs(globIdx_receiveDirichlet) = 0;
              obj.ls.lhs(globIdx_receiveDirichlet, :) = 0.0;
              obj.ls.lhs(globIdx_receiveDirichlet, globIdx_receiveDirichlet) = -1;
              obj.ls.lhs(globIdx_receiveDirichlet, globIdx_sendDirichlet) = coeffs;
          end
          % obj.problemHS receives NEUMANN from Part
          for inode=hsGammaNodes
              posNode = obj.problemHS.mesh.pos( inode );
              posNode_other = obj.problemHS.mesh.posFixed(inode) - obj.problemPart.mesh.shiftFixed;
              % Compute normal
              % Find parentEl
              ownerEls = obj.meshSubdomain.conPointCell.getLocalCon( inode );
              idxParentEl = ownerEls( logical(obj.meshSubdomain.activeElements(ownerEls)) );
              parentEl = obj.meshSubdomain.getElement( idxParentEl );
              normal = posNode - parentEl.centroid;
              normal = normal / abs(normal);
              % Get corresponding element neighbour mesh
              iowner = obj.problemPart.mesh.findOwnerElement( posNode_other );
              e = obj.problemPart.mesh.getElement( iowner );
              %TODO: THink about conductivities here!
              dn_coeffs = normal*[e.gradBaseFuns{1}(posNode_other), e.gradBaseFuns{2}(posNode_other)];

              globIdx_receiveNeumann = obj.problemHS.numbering( inode );
              globIdx_sendNeumann = obj.problemPart.numbering( e.con );
              
              obj.ls.lhs(globIdx_receiveNeumann, globIdx_sendNeumann) = obj.ls.lhs(globIdx_receiveNeumann, globIdx_sendNeumann) - ...
                  obj.problemHS.k * dn_coeffs;
          end

          obj.problemPart.assembleInactiveNodes(obj.ls);
          obj.problemHS.assembleInactiveNodes(obj.ls);
          %% Solve
          obj.ls.solve();
          %% post-iteration
          obj.problemPart.U = obj.ls.U(obj.problemPart.numbering);
          obj.problemHS.U = obj.ls.U(obj.problemHS.numbering);
          % Interpolate inactive nodes of obj.problemPart
          obj.problemPart.interpolateInactive( obj.problemHS );
          obj.problemHS.interpolateInactive( obj.problemPart );
          % Save interface positions for post
          obj.posInterface = obj.problemPart.mesh.posFixed( partGammaNodes );
        end
        function [t] = getTime(obj)
            t = obj.problemHS.time;
        end
    end
end

function plotActiveInactive( scheme )
    arguments
        scheme MyDDScheme
    end
    figure
    scatter(scheme.meshPart.posFixed(~logical(scheme.meshPart.activeNodes)), 1, "r")
    hold on
    scatter(scheme.meshPart.posFixed(logical(scheme.meshPart.activeNodes)), 1, "k")

    scatter(scheme.meshSubdomain.posFixed(~logical(scheme.meshSubdomain.activeNodes)), 1.25, "r")
    scatter(scheme.meshSubdomain.posFixed(logical(scheme.meshSubdomain.activeNodes)), 1.25, "k")
end