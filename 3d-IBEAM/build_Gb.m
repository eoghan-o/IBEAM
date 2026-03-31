function Gb_T = build_Gb(local_face, n_faces, vertices, a, b)
    n_v  = size(vertices, 1);
    Gb_T = zeros(n_faces, n_v);

    for e = 1:n_faces
        i     = local_face(e, 1);
        j     = local_face(e, end);
        e_vec = (vertices(j,:) - vertices(i,:))';  % edge vector v_j - v_i, like 2D
        p_e   = a^(-1) * b' * e_vec;
        Gb_T(e, i) = -bernoulli( p_e);
        Gb_T(e, j) =  bernoulli(-p_e);
    end
end