function [local_edge, local_face] = build_edge_and_face_list(vertices)
    % Returns:
    % local_edge - (n_edges x 2): [i j] pairs for G and Gb (real + fake edges)
    % local_face - (n_faces x 3): triangular faces for W (convex hull)

    n_v = size(vertices, 1);

    % local_face: convex hull triangulation for W
    K          = convhull(vertices(:,1), vertices(:,2), vertices(:,3));
    local_face = K;  % (n_real x 3)
    n_real     = size(local_face, 1);

    % Real edges: one per face, using first and last vertex (like before)
    real_edges = [local_face(:,1), local_face(:,end)];

    % Fake edges: all vertex pairs not already a real edge
    all_pairs  = nchoosek(1:n_v, 2);
    n_all      = size(all_pairs, 1);
    real_set   = sort(real_edges, 2);
    fake_edges = zeros(n_all, 2);
    n_fake     = 0;
    for i = 1:n_all
        pair = sort(all_pairs(i,:));
        if ~ismember(pair, real_set, 'rows')
            n_fake = n_fake + 1;
            fake_edges(n_fake,:) = all_pairs(i,:);
        end
    end
    fake_edges = fake_edges(1:n_fake,:);
    local_edge = [real_edges; fake_edges];
end