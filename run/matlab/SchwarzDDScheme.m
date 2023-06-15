classdef SchwarzDDScheme < handle
    properties
        ls LinearSystem
        numberingPart
        numberingSubdomain
        meshPart Mesh
        meshSubdomain Mesh
        problemPart Subproblem
        problemHS Subproblem
        posInterface = [];
        pad = 0.0;
    end
    methods
        function [obj] = SchwarzDDScheme(params)
            %MESHING
            nelsPart = (params.rightBound - params.leftBound) * params.meshDensity;
            nelsHS = params.radiusSubdomain*2*params.meshDensity;
            obj.meshPart = Mesh( params.leftBound, params.rightBound, nelsPart );
            obj.meshSubdomain = Mesh( params.x0-params.radiusSubdomain, params.x0+params.radiusSubdomain, nelsHS );
            obj.pad = params.pad;
            
            ndofs = obj.meshPart.nnodes + obj.meshSubdomain.nnodes;
            %NUMBERING DOFS
            obj.numberingPart = 1:obj.meshPart.nnodes;
            obj.numberingSubdomain = (obj.numberingPart(end)+1):(obj.numberingPart(end)+obj.meshSubdomain.nnodes);
            obj.ls = LinearSystem( ndofs );
            
            % Subproblems
            params.vadv = -params.speed;
            params.vsource = 0.0;
            params.vdomain = params.speed;
            obj.problemHS = Subproblem(obj.meshSubdomain, obj.numberingSubdomain, params);
            params.vadv = 0.0;
            params.vsource = params.speed;
            params.vdomain = 0.0;
            params.bcNodes = [1, obj.meshPart.nnodes];
            obj.problemPart = Subproblem(obj.meshPart, obj.numberingPart, params);
 
        end
        function [obj] = iterate(obj)
          %% pre-iteration
          obj.ls.cleanup();
          % Domain OPS
          obj.problemHS.mesh.intersect( obj.problemPart.mesh, true );
          obj.problemPart.preIterate();
          obj.problemHS.preIterate();
          obj.problemHS.mesh.intersect( obj.problemPart.mesh, false );
          obj.problemPart.mesh.substract( obj.problemHS.mesh, true, obj.pad );
          % Determine Gamma boundary
          obj.problemHS.updateGammaNodes( obj.meshPart, true );
          obj.problemPart.updateGammaNodes( obj.meshSubdomain, true );

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
          % obj.problemHS receives Dirichlet from Part
          for inode=obj.problemHS.gammaNodes
              posNode_other = obj.problemHS.mesh.posFixed(inode) - obj.problemPart.mesh.shiftFixed;
              iowner = obj.problemPart.mesh.findOwnerElement( posNode_other );
              e = obj.problemPart.mesh.getElement( iowner );
              coeffs = [e.baseFuns{1}(posNode_other), e.baseFuns{2}(posNode_other)];

              globIdx_receiveDirichlet = obj.problemHS.numbering( inode );
              globIdx_sendDirichlet = obj.problemPart.numbering( e.con );

              obj.ls.rhs(globIdx_receiveDirichlet) = 0;
              obj.ls.lhs(globIdx_receiveDirichlet, :) = 0.0;
              obj.ls.lhs(globIdx_receiveDirichlet, globIdx_receiveDirichlet) = -1;
              obj.ls.lhs(globIdx_receiveDirichlet, globIdx_sendDirichlet) = coeffs;
          end

          obj.problemPart.assembleInactiveNodes(obj.ls);
          obj.problemHS.assembleInactiveNodes(obj.ls);
          %% Solve
          obj.ls.solve();
          %% post-iteration
          obj.problemPart.U = obj.ls.U(obj.problemPart.numbering);
          obj.problemHS.U = obj.ls.U(obj.problemHS.numbering);
          % Interpolate inactive nodes of obj.problemPart
          obj.problemPart.interpolateOther( obj.problemHS );
          obj.problemHS.interpolateInactive( obj.problemPart );
          % Save interface positions for post
          obj.posInterface = obj.problemPart.mesh.posFixed( obj.problemPart.gammaNodes );
        end
        function [t] = getTime(obj)
            t = obj.problemHS.time;
        end
        function [] = plotInterface(obj)
            posGammaPart = obj.problemPart.mesh.pos( obj.problemPart.gammaNodes );
            posGammaHS = obj.problemHS.mesh.pos( obj.problemHS.gammaNodes ) + obj.problemHS.mesh.shiftFixed;
            for idx=1:length(posGammaPart)
                xGamma = posGammaPart(idx);
                if idx==1
                    xline(xGamma, ...
                '--r', 'DisplayName', "Schwarz, $\Gamma$ fixed domain")
                else
                    xline(xGamma, ...
                '--r', 'HandleVisibility', "off")
                end
            end
            for idx=1:length(posGammaHS)
                xGamma = posGammaHS(idx);
                if idx==1
                    xline(xGamma, ...
                '--k', 'DisplayName', "Schwarz, $\Gamma$ moving subdomain")
                else
                    xline(xGamma, ...
                '--k', 'HandleVisibility', "off")
                end
            end
        end
    end
end

function plotActiveInactive( scheme )
    arguments
        scheme SchwarzDDScheme
    end
    figure
    scatter(scheme.meshPart.posFixed(~logical(scheme.meshPart.activeNodes)), 1, "r")
    hold on
    scatter(scheme.meshPart.posFixed(logical(scheme.meshPart.activeNodes)), 1, "k")

    scatter(scheme.meshSubdomain.posFixed(~logical(scheme.meshSubdomain.activeNodes)), 1.25, "r")
    scatter(scheme.meshSubdomain.posFixed(logical(scheme.meshSubdomain.activeNodes)), 1.25, "k")
end