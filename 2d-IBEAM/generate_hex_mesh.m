function [elements, nodes, n_nodes] = generate_hex_mesh(n)
% Regular hexagonal mesh on [0,1]^2.
% n controls number of hexagons across the width.

    % Hex geometry: pointy-top orientation
    % This gives cleaner boundary alignment on [0,1]^2
    r  = 1 / (2 * (n - 1));          % circumradius
    dx = r * sqrt(3);                  % horizontal center spacing
    dy = 2 * r;                        % vertical center spacing (pointy-top)

    snap_tol = r * 0.4;
    tol      = 1e-10;

    node_map = containers.Map('KeyType','char','ValueType','double');
    nodes    = zeros(n^2 * 6, 2);
    n_nodes  = 0;

    function idx = get_or_add_node(x, y)
        if x < snap_tol,           x = 0; elseif x > 1-snap_tol, x = 1; end
        if y < snap_tol,           y = 0; elseif y > 1-snap_tol, y = 1; end
        key = sprintf('%.8f_%.8f', round(x/tol)*tol, round(y/tol)*tol);
        if isKey(node_map, key)
            idx = node_map(key);
        else
            n_nodes          = n_nodes + 1;
            nodes(n_nodes,:) = [x, y];
            node_map(key)    = n_nodes;
            idx              = n_nodes;
        end
    end

    % Pointy-top hex: CCW corners at 90, 150, 210, 270, 330, 30 degrees
    angles = pi/2 + (0:5) * pi/3;

    n_cols = round(1 / dx) + 2;
    n_rows = round(1 / (dy * 0.75)) + 2;

    elements = cell(n_cols * n_rows, 1);
    e = 1;

    for j = 0:n_cols
        cx     = j * dx + mod(0,2);           % no column offset needed
        for i = 0:n_rows
            % Row stagger: odd rows offset by dx/2
            cy = i * dy * 0.75;
            if mod(i, 2) == 1
                cxi = cx + dx/2;
            else
                cxi = cx;
            end

            % Skip centers too far outside domain
            if cxi < -dx || cxi > 1+dx, continue; end
            if cy  < -dy || cy  > 1+dy, continue; end

            conn = zeros(1,6);
            for k = 1:6
                conn(k) = get_or_add_node(cxi + r*cos(angles(k)), ...
                                          cy  + r*sin(angles(k)));
            end

            [~, keep_idx] = unique(conn, 'stable');
            conn = conn(sort(keep_idx));

            if length(conn) >= 3
                elements{e} = conn;
                e = e + 1;
            end
        end
    end

    elements = elements(1:e-1);
    nodes    = nodes(1:n_nodes,:);
    n_nodes  = size(nodes, 1);
end