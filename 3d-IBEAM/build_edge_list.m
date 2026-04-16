function local_edge = build_edge_list(vertices)

    n_v = size(vertices, 1);

    % ------------------------------------------------------------
    % STEP 1: Define REAL edges purely geometrically
    % (use convex hull triangles, but safely)
    % ------------------------------------------------------------
    K = convhulln(vertices);  % <-- IMPORTANT: use convhulln (true 3D)

    all_edges = [K(:,1) K(:,2);
                 K(:,2) K(:,3);
                 K(:,3) K(:,1)];

    real_edges = unique(sort(all_edges, 2), 'rows');

    % ------------------------------------------------------------
    % STEP 2: Add all missing pairs (complete graph closure)
    % ------------------------------------------------------------
    all_pairs = nchoosek(1:n_v, 2);

    is_real = ismember(all_pairs, real_edges, 'rows');
    fake_edges = all_pairs(~is_real, :);

    % ------------------------------------------------------------
    % STEP 3: Combine
    % ------------------------------------------------------------
    local_edge = [real_edges; fake_edges];

end