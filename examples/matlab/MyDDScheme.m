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
        currRadiusSubdomain = -1;
        hasChangedRadiusSubdomain = false;
        tstepcounter = 0.0;
        label = "My scheme";
    end
    methods
        function [obj] = MyDDScheme(params)
            %MESHING
            obj.currRadiusSubdomain = params.radiusSubdomain;
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
            params.bcNodes = [1, obj.meshPart.nnodes];
            obj.problemPart = Subproblem(obj.meshPart, obj.numberingPart, params);
            params.vadv = -params.speed;
            params.vsource = 0.0;
            params.vdomain = params.speed;
            params.bcNodes = [];
            obj.problemHS = Subproblem(obj.meshSubdomain, obj.numberingSubdomain, params);
        end
        function [obj] = iterate(obj)
          obj.tstepcounter = obj.tstepcounter + 1;
          %% pre-iteration
          obj.ls.cleanup();
          % Motion is done here
          obj.problemHS.mesh.intersect( obj.problemPart.mesh, true );
          obj.problemPart.preIterate();
          obj.problemHS.preIterate();
          % Domain OPS
          if (obj.hasChangedRadiusSubdomain)
            obj.problemHS.mesh.intersectBall( obj.problemHS.x0, obj.currRadiusSubdomain, false );
          end
          obj.problemHS.mesh.intersect( obj.problemPart.mesh, false );
          obj.problemPart.mesh.substract( obj.problemHS.mesh );
          % Determine Gamma boundary
          % Determine Gamma boundary
          obj.problemHS.updateGammaNodes( obj.meshPart, false );
          obj.problemPart.updateGammaNodes( obj.meshSubdomain, false );

          %% Assembly
          obj.problemPart.assemblePDE(obj.ls);
          obj.problemHS.assemblePDE(obj.ls);
          % INTERFACE
          % obj.problemPart receives DIRICHLET from HS
          for inode=obj.problemPart.gammaNodes
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
          for inode=obj.problemHS.gammaNodes
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
              
              obj.ls.lhs(globIdx_receiveNeumann, globIdx_sendNeumann) = obj.ls.lhs(globIdx_receiveNeumann, globIdx_sendNeumann) ...
                  - obj.problemHS.k * dn_coeffs;
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
        end
        function [t] = getTime(obj)
            t = obj.problemHS.time;
        end
        function [obj] = setDt(obj, dt)
            obj.problemHS.dt = dt;
            obj.problemPart.dt = dt;
        end
        function [] = plotInterface(obj)
            posGammaPart = obj.problemPart.mesh.pos( obj.problemPart.gammaNodes );
            for idx=1:length(obj.problemPart.gammaNodes)
                xGamma = posGammaPart(idx);
                if idx==1
                    xline(xGamma, ...
                'DisplayName', obj.label + ", $\Gamma$")
                else
                    xline(xGamma, ...
                'HandleVisibility', "off")
                end
            end
        end
        function obj = setRadiusSubdomain(obj, r)
            obj.currRadiusSubdomain = r;
            obj.hasChangedRadiusSubdomain = true;
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