classdef Mesh  < handle
    % 1D mesh
  properties
    nels
    nnodes
    pos
    shiftFixed = 0.0 %pos + shiftFixed = posFixed
    posFixed %position in fixed reference frame
    conCellPoint Connectivity
    conPointCell Connectivity
    bounNodes %indices of boundary nodes
    mass
    left = 1;
    right = -1;
    activeElements
    activeNodes
  end
  methods
      function obj = Mesh( left, right, nels )
          % MESH Builds 1D mesh
          %    mesh = MESH( left, right, nels )
          obj.left = left;
          obj.right = right;
          obj.nels = nels;
          obj.nnodes = nels + 1;
          obj.pos = linspace( left, right, obj.nnodes )';
          obj.posFixed = obj.pos;
          obj.conCellPoint = Connectivity( 1, 0, obj.nels, obj.nnodes, [(1:obj.nels)', (2:(obj.nels+1))']);
          obj.conPointCell = obj.conCellPoint.transpose();
          obj.activeElements = ones([obj.nels, 1]);
          obj.updateActiveNodes();
          obj.postActivation();
      end
      function [obj, newInterface] = substract( obj, otherMesh, reset, pad )
          % 1D
          % If element is contained in other mesh, deactivate
          arguments
              obj Mesh
              otherMesh Mesh
              reset = true
              pad = 0.0;
          end
          if reset
            obj.activeElements = ones([obj.nels, 1]);
          end
          % Assuming otherMesh has no holes
          otherMeshLeft = otherMesh.posFixed( find(otherMesh.activeNodes, 1) );
          otherMeshLeft = otherMeshLeft + pad;
          otherMeshRight = otherMesh.posFixed( find(otherMesh.activeNodes, 1, "last") );
          otherMeshRight = otherMeshRight - pad;
          for iel=1:obj.nels
              e = obj.getElement(iel);
              elLeft = e.pos(1);
              elRight = e.pos(2);
              
              if (otherMeshLeft <= elLeft ) && (otherMeshRight >= elRight)
                  obj.activeElements(iel) = 0.0;
              end
          end
          [~, newInterface] = obj.postActivateElements();
      end
      function [obj, newInterface] = intersect( obj, otherMesh, reset )
          % 1D
          % If element is NOT contained in other mesh, deactivate
          arguments
              obj Mesh
              otherMesh Mesh
              reset = true
          end
          if reset
            obj.activeElements = ones([obj.nels, 1]);
          end
          for iel=1:obj.nels
              e = obj.getElement(iel);
              elLeft = e.pos(1) + obj.shiftFixed;
              elRight = e.pos(2) + obj.shiftFixed;

              if (otherMesh.posFixed(1) > elLeft ) || (otherMesh.posFixed(end) < elRight)
                  obj.activeElements(iel) = 0;
              end
          end
          [~, newInterface] = obj.postActivateElements();
      end
      function [obj, newInterface] = intersectBall( obj, center, radius, reset )
          % 1D
          % If element is NOT contained in other mesh, deactivate
          arguments
              obj Mesh
              center
              radius
              reset = true
          end
          if reset
            obj.activeElements = ones([obj.nels, 1]);
          end
          leftBound = center - radius;
          rightBound = center + radius;
          for iel=1:obj.nels
              e = obj.getElement(iel);
              elLeft = e.pos(1);
              elRight = e.pos(2);

              if (leftBound > elLeft ) || (rightBound < elRight)
                  obj.activeElements(iel) = 0;
              end
          end
          [~, newInterface] = obj.postActivateElements();
      end
      function [obj, newInterface] = postActivateElements(obj)
          olbBoundary = obj.bounNodes;
          obj.updateActiveNodes();
          obj.postActivation();
          newInterface = obj.bounNodes( find( ~ismember(obj.bounNodes, olbBoundary) ) );
      end
      function obj = postActivation(obj)
        % POS-ACTIVATION
        obj.computeMass();%RECOMPUTE MASS MATRIX
        obj.findBoundary();
      end
      function obj = updateActiveNodes(obj)
          obj.activeNodes = zeros([obj.nnodes, 1]);
          for iel=1:obj.nels
              if obj.activeElements(iel)
                  obj.activeNodes( obj.conCellPoint.getLocalCon(iel) ) = 1;
              end
          end
      end
      function obj = findBoundary(obj)
          obj.bounNodes = [];
          for inode=1:obj.nnodes
              locCon = obj.conPointCell.getLocalCon( inode );
              numIncidentActiveEls = sum( obj.activeElements(locCon));
              if (numIncidentActiveEls == 1)
                  obj.bounNodes = [obj.bounNodes, inode];
              end
          end
      end
      function obj = computeMass(obj)
          i = [];
          j = [];
          v = [];
          for ielem=1:obj.nels
              if ~(obj.activeElements(ielem))
                  continue;
              end
              % Load element
              e = obj.getElement( ielem ); 
              % Compute contributions
              Mloc = zeros([e.nnodes, e.nnodes]);
              for igp=1:e.ngpoints
                  for inode=1:e.nnodes
                      for jnode=1:e.nnodes
                          Mloc(inode, jnode) = Mloc(inode, jnode) + ...
                              e.baseFuns{inode}(e.gpos(igp))*e.baseFuns{jnode}(e.gpos(igp))*...
                              e.gpweight(igp)*e.vol;
                      end
                  end
              end
              for inode=1:e.nnodes
                  for jnode=1:e.nnodes
                      i = [i, e.con(inode)];
                      j = [j, e.con(jnode)];
                      v = [v, Mloc(inode, jnode)];
                  end
              end
          end
          obj.mass = sparse( i, j, v, obj.nnodes, obj.nnodes );
      end
      function element = getElement(obj, ielem)
          locCon = obj.conCellPoint.getLocalCon( ielem );
          element = LineElement(obj.pos( locCon ), locCon );
      end
      function [obj] = updatePosFixed( obj, shift )
          obj.shiftFixed = shift;
          obj.posFixed = obj.pos + obj.shiftFixed;
      end
      function idxOwnerEl = findOwnerElement(obj, x)
          idxOwnerEl = -1;
          for ielem=1:obj.nels
              if (~obj.activeElements(ielem))
                  continue;
              end
              e = obj.getElement(ielem);
              locPos = obj.pos( e.con );
              if (locPos(1) <= x) && (locPos(2) >= x)
                  idxOwnerEl = ielem;
                  break;
              end
          end
      end
  end
end
