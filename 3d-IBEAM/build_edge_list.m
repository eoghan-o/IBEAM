function local_edge = build_edge_list(vertices)
    n_v = size(vertices, 1);

    % Real edges from convex hull faces
    K          = convhull(vertices(:,1), vertices(:,2), vertices(:,3));
    % Extract unique edges from triangular faces
    all_edges  = [K(:,1) K(:,2); K(:,2) K(:,3); K(:,1) K(:,3)];
    real_edges = unique(sort(all_edges, 2), 'rows');

    % Fake edges: all vertex pairs not already real edges
    all_pairs  = nchoosek(1:n_v, 2);
    fake_edges = zeros(size(all_pairs,1), 2);
    n_fake     = 0;
    for i = 1:size(all_pairs,1)
        if ~ismember(all_pairs(i,:), real_edges, 'rows')
            n_fake = n_fake + 1;
            fake_edges(n_fake,:) = all_pairs(i,:);
        end
    end
    local_edge = [real_edges; fake_edges(1:n_fake,:)];
end