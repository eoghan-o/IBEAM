function local_face = build_face_list(vertices)
    if size(vertices,1) == 8
    local_face = [
            1 2 3 4;
            5 6 7 8;
            1 2 6 5;
            2 3 7 6;
            3 4 8 7;
            4 1 5 8
        ];
    return;
    end

    % Compute convex hull to get real faces
    % Each row of real_faces is a triangle of vertex indices (local)
    try
        K = convhull(vertices(:,1), vertices(:,2), vertices(:,3));
        real_faces = K;  % (n_real x 3) triangular faces
    catch
        error('convhull failed for element vertices');
    end
    n_v      = size(vertices, 1);
    n_real   = size(real_faces, 1);
    
    % Virtual faces: all triples not already in real_faces
    % (analogous to non-adjacent edges in 2D)
    all_triples = nchoosek(1:n_v, 3);
    n_all    = size(all_triples, 1);
    fake_faces = zeros(n_all, 3);
    n_fake   = 0;
    real_set = sort(real_faces, 2);
    for i = 1:n_all
        triple = sort(all_triples(i,:));
        if ~ismember(triple, real_set, 'rows')
            n_fake = n_fake + 1;
            fake_faces(n_fake,:) = all_triples(i,:);
        end
    end
    fake_faces = fake_faces(1:n_fake, :);
    local_face = [real_faces; fake_faces];

end