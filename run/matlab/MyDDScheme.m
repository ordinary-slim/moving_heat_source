classdef MyDDScheme
    properties
        ls LinearSystem
        numberingPart
        numberingSubdomain
        meshPart Mesh
        meshSubdomain Mesh
        problemPart Subproblem
        problemHS Subproblem
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
          obj.problemPart.preIterate();
          obj.problemHS.preIterate();
          % Domain OPS
          [~, partGammaNodes] = obj.problemPart.mesh.substract( obj.problemHS.mesh );
          obj.problemHS.mesh.intersect( obj.problemPart.mesh );
          hsGammaNodes = obj.problemHS.mesh.bounNodes;%TODO: fix this
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
              parentEl = obj.meshSubdomain.getElement( obj.meshSubdomain.conPointCell.getLocalCon( inode ) );
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
          inactiveNodesPart = find(~obj.problemPart.mesh.activeNodes);
          for inode=inactiveNodesPart'
              posNode_other = obj.problemPart.mesh.posFixed(inode) - obj.problemHS.mesh.shiftFixed;
              obj.problemPart.U(inode) =  obj.problemHS.evaluate( posNode_other );
          end
        end
    end
end
