fprintf('=== Gb vs G diagnostic ===\n\n');

for k = 2:4
    h = 1/k;
    [nodes, elements] = generate_mesh3d(h);

    T        = 1;
    indices  = elements{T};
    vertices = nodes(indices,:);
    n_v      = size(vertices,1);

    local_face = build_face_list(vertices);
    n_faces    = size(local_face,1);

    G_T  = build_G(local_face, n_v, n_faces);
    Gb_T = build_Gb(local_face, n_faces, vertices, 1, [0;0;0]);

    fprintf('k=%d  norm(Gb-G) = %.2e\n', k, norm(Gb_T - G_T, 'fro'));

    % Now check: in 2D, what is p_e for edge (i,j)?
    % p_e = a^{-1} * b' * e_vec  where e_vec = v_j - v_i
    % In 3D we use face normal instead - is that correct?
    % Let's print p values for each face
    fprintf('  Face  p_face  bernoulli(p)  bernoulli(-p)\n');
    b_test = [0;0;1];
    for e = 1:n_faces
        i     = local_face(e,1);
        j     = local_face(e,end);
        vi    = vertices(i,:);
        vj    = vertices(j,:);
        % 3D current: face normal
        e1     = vertices(local_face(e,2),:) - vi;
        e2     = vertices(local_face(e,3),:) - vi;
        normal = cross(e1, e2);
        len    = norm(normal);
        normal_unit = normal/len;
        p_face = b_test' * normal_unit' * len;
        % What it should be: edge vector v_j - v_i like in 2D
        e_vec  = (vj - vi)';
        p_edge = b_test' * e_vec;
        fprintf('  e=%d  p_face=%.4f  p_edge=%.4f\n', e, p_face, p_edge);
    end
    fprintf('\n');
end