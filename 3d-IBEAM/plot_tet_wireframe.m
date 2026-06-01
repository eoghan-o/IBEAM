function plot_tet_wireframe(nodes, elements, varargin)
% PLOT_TET_WIREFRAME  Render a tetrahedral mesh as a wireframe
%
% Usage:
%   plot_tet_wireframe(nodes, elements)
%   plot_tet_wireframe(nodes, elements, 'overlay', true)
%   plot_tet_wireframe(nodes, elements, 'color', [0.1 0.1 0.6], 'linewidth', 0.8)
%
% Inputs:
%   nodes    - (n_nodes x 3) node coordinates
%   elements - cell array of (4 x 1) vertex index lists (your format)

    % --- Parse optional args ---
    p = inputParser;
    addParameter(p, 'overlay',   false,          @islogical);
    addParameter(p, 'color',     [0.15 0.15 0.15], @isnumeric);
    addParameter(p, 'facecolor', [0.82 0.87 0.94], @isnumeric);
    addParameter(p, 'linewidth', 0.6,            @isnumeric);
    addParameter(p, 'alpha',     0.15,           @isnumeric);
    parse(p, varargin{:});
    opt = p.Results;

    % --- Build unique edge list from tetrahedra ---
    % Each tet has 6 edges: all pairs from 4 nodes
    tet_edge_pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    n_tets = numel(elements);
    all_edges = zeros(n_tets * 6, 2);

    for t = 1:n_tets
        verts = elements{t};   % 4 node indices
        for e = 1:6
            i1 = verts(tet_edge_pairs(e,1));
            i2 = verts(tet_edge_pairs(e,2));
            all_edges((t-1)*6 + e, :) = sort([i1, i2]);
        end
    end

    % Remove duplicate edges (interior edges shared by multiple tets)
    unique_edges = unique(all_edges, 'rows');

    % --- Plot ---
    figure;
    hold on;

    if opt.overlay
        % Build face list for solid surface (4 faces per tet)
        tet_face_triples = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
        all_faces = zeros(n_tets * 4, 3);
        for t = 1:n_tets
            verts = elements{t};
            for f = 1:4
                all_faces((t-1)*4 + f, :) = verts(tet_face_triples(f,:));
            end
        end
        % Solid mesh underneath
        patch('Vertices', nodes, 'Faces', all_faces, ...
              'FaceColor',  opt.facecolor, ...
              'FaceAlpha',  opt.alpha, ...
              'EdgeColor',  'none', ...
              'FaceLighting', 'gouraud', ...
              'AmbientStrength', 0.5);
        camlight; lighting gouraud;
    end

    % --- Draw all edges in one vectorized call (fast) ---
    % Build NaN-separated line segments: [x1 x2 NaN; ...] etc.
    e1 = unique_edges(:,1);
    e2 = unique_edges(:,2);
    nan_col = nan(size(e1));

    Xw = [nodes(e1,1), nodes(e2,1), nan_col]';   % 3 x n_edges
    Yw = [nodes(e1,2), nodes(e2,2), nan_col]';
    Zw = [nodes(e1,3), nodes(e2,3), nan_col]';

    plot3(Xw(:), Yw(:), Zw(:), ...
          'Color',     opt.color, ...
          'LineWidth', opt.linewidth);

    axis equal; grid on; view(3);
    xlabel('x'); ylabel('y'); zlabel('z');
    title(sprintf('Tet mesh — %d nodes, %d elements, %d edges', ...
          size(nodes,1), n_tets, size(unique_edges,1)));
    hold off;
end