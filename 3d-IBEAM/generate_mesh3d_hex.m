function [nodes, elements] = generate_mesh3d_hex(h)
    % Generate a structured hexahedral mesh on [0,1]^3
    % Each cube cell is a single 8-node element
    %
    % Input:
    %   h        - mesh spacing (same in x, y, z)
    % Outputs:
    %   nodes    - (n_nodes x 3) node coordinates
    %   elements - cell array of (8 x 1) vertex index lists

    n  = round(1/h);
    x  = linspace(0, 1, n+1);
    y  = linspace(0, 1, n+1);
    z  = linspace(0, 1, n+1);

    % Node indexing: node_id(i,j,k) = i + j*(n+1) + k*(n+1)^2
    node_id = @(i,j,k) i + j*(n+1) + k*(n+1)^2;

    % Build node list
    n_nodes = (n+1)^3;
    nodes   = zeros(n_nodes, 3);
    for k = 0:n
        for j = 0:n
            for i = 0:n
                id          = node_id(i, j, k) + 1;  % 1-indexed
                nodes(id,:) = [x(i+1), y(j+1), z(k+1)];
            end
        end
    end

    % Vertex ordering (same as your tetrahedral version):
    % v1=(i,j,k)   v2=(i+1,j,k)   v3=(i+1,j+1,k)  v4=(i,j+1,k)
    % v5=(i,j,k+1) v6=(i+1,j,k+1) v7=(i+1,j+1,k+1) v8=(i,j+1,k+1)

    n_cubes  = n^3;
    elements = cell(n_cubes, 1);
    el_count = 0;

    for k = 0:n-1
        for j = 0:n-1
            for i = 0:n-1
                % 8 corners of this cube (1-indexed)
                v(1) = node_id(i,   j,   k  ) + 1;
                v(2) = node_id(i+1, j,   k  ) + 1;
                v(3) = node_id(i+1, j+1, k  ) + 1;
                v(4) = node_id(i,   j+1, k  ) + 1;
                v(5) = node_id(i,   j,   k+1) + 1;
                v(6) = node_id(i+1, j,   k+1) + 1;
                v(7) = node_id(i+1, j+1, k+1) + 1;
                v(8) = node_id(i,   j+1, k+1) + 1;

                el_count = el_count + 1;
                elements{el_count} = v(:);  % store as column vector
            end
        end
    end
end