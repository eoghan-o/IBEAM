function total = polyintegrate3d(v, f, n_quad, integrand)
    dummy = integrand(v(1,:));
    total = zeros(size(dummy));

    anchor = v(1,:);

    for i = 1:length(f)
        face = f{i};
        n_f = length(face);

        if any(face == 1)
            continue
        end

        for j = 1:(n_f - 2)
            a = v(face(1),   :);
            b = v(face(j+1), :);
            c = v(face(j+2), :);

            tet = [anchor; a; b; c];

            [xq, wq] = tetquad(tet, n_quad);

            n_pts = size(xq, 1);
            for p = 1:n_pts
                total = total + wq(p) * integrand(xq(p,:));
            end
        end
    end
end