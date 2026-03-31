function [nodes, elements] = generate_mesh3d(h)
    % Generate a structured tetrahedral mesh on [0,1]^3
    % Each cube cell is split into 6 tetrahedra
    %
    % Input:
    %   h        - mesh spacing (same in x, y, z)
    % Outputs:
    %   nodes    - (n_nodes x 3) node coordinates
    %   elements - cell array of (4 x 1) vertex index lists

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
                id         = node_id(i, j, k) + 1;  % 1-indexed
                nodes(id,:) = [x(i+1), y(j+1), z(k+1)];
            end
        end
    end

    % Build element list: split each cube into 6 tets
    % Standard 6-tet decomposition of a unit cube with corners
    % v1..v8 ordered as:
    %   v1=(i,j,k)   v2=(i+1,j,k)   v3=(i+1,j+1,k)  v4=(i,j+1,k)
    %   v5=(i,j,k+1) v6=(i+1,j,k+1) v7=(i+1,j+1,k+1) v8=(i,j+1,k+1)
    %
    % The 6 tets use the Kuhn/Freudenthal decomposition which is
    % consistent across neighbouring cubes (no hanging nodes):
    tet_local = [1 2 4 5;
                 2 4 5 8;
                 2 5 6 8;
                 2 3 4 8;
                 2 6 7 8;
                 2 3 7 8];

    n_cubes   = n^3;
    elements  = cell(6 * n_cubes, 1);
    el_count  = 0;

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

                for t = 1:6
                    el_count = el_count + 1;
                    elements{el_count} = v(tet_local(t,:));
                end
            end
        end
    end
end