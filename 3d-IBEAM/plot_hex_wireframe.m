function plot_hex_wireframe(nodes, elements, varargin)
% PLOT_HEX_WIREFRAME  Render a hexahedral mesh as a wireframe
%
% Usage:
%   plot_hex_wireframe(nodes, elements)
%   plot_hex_wireframe(nodes, elements, 'overlay', true)
%   plot_hex_wireframe(nodes, elements, 'color', [0.1 0.1 0.6], 'linewidth', 0.8)
%
% Inputs:
%   nodes    - (n_nodes x 3) node coordinates
%   elements - cell array of (8 x 1) vertex index lists (your format)

    % --- Parse optional args ---
    p = inputParser;
    addParameter(p, 'overlay',   false,            @islogical);
    addParameter(p, 'color',     [0.15 0.15 0.15], @isnumeric);
    addParameter(p, 'facecolor', [0.82 0.87 0.94], @isnumeric);
    addParameter(p, 'linewidth', 0.6,              @isnumeric);
    addParameter(p, 'alpha',     0.15,             @isnumeric);
    parse(p, varargin{:});
    opt = p.Results;

    % --- Hex topology ---
    % 12 edges per hex (pairs of local vertex indices)
    hex_edge_pairs = [1 2; 2 3; 3 4; 4 1; ...   % bottom face
                      5 6; 6 7; 7 8; 8 5; ...   % top face
                      1 5; 2 6; 3 7; 4 8];      % vertical edges

    % 6 faces per hex (quads, used for overlay)
    hex_face_quads = [1 2 3 4; ...   % bottom
                      5 6 7 8; ...   % top
                      1 2 6 5; ...   % front
                      3 4 8 7; ...   % back
                      2 3 7 6; ...   % right
                      1 4 8 5];     % left

    % --- Build unique edge list ---
    n_els     = numel(elements);
    all_edges = zeros(n_els * 12, 2);

    for el = 1:n_els
        verts = elements{el};
        for e = 1:12
            i1 = verts(hex_edge_pairs(e,1));
            i2 = verts(hex_edge_pairs(e,2));
            all_edges((el-1)*12 + e, :) = sort([i1, i2]);
        end
    end

    unique_edges = unique(all_edges, 'rows');

    % --- Plot ---
    figure;
    hold on;

    if opt.overlay
        % Build quad face list for solid surface
        all_faces = zeros(n_els * 6, 4);
        for el = 1:n_els
            verts = elements{el};
            for f = 1:6
                all_faces((el-1)*6 + f, :) = verts(hex_face_quads(f,:))';
            end
        end
        patch('Vertices', nodes, 'Faces', all_faces, ...
              'FaceColor',      opt.facecolor, ...
              'FaceAlpha',      opt.alpha, ...
              'EdgeColor',      'none', ...
              'FaceLighting',   'gouraud', ...
              'AmbientStrength', 0.5);
        camlight; lighting gouraud;
    end

    % --- Vectorized plot3 with NaN separators ---
    e1 = unique_edges(:,1);
    e2 = unique_edges(:,2);
    nan_col = nan(size(e1));

    Xw = [nodes(e1,1), nodes(e2,1), nan_col]';
    Yw = [nodes(e1,2), nodes(e2,2), nan_col]';
    Zw = [nodes(e1,3), nodes(e2,3), nan_col]';

    plot3(Xw(:), Yw(:), Zw(:), ...
          'Color',     opt.color, ...
          'LineWidth', opt.linewidth);

    axis equal; grid on; view(3);
    xlabel('x'); ylabel('y'); zlabel('z');
    title(sprintf('Hex mesh — %d nodes, %d elements, %d edges', ...
          size(nodes,1), n_els, size(unique_edges,1)));
    hold off;
end