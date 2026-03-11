function W = assemble_weight_matrix(vertices)
    w = zeros(6,1);

    v1 = vertices([1 2 3], :);
    v2 = vertices([1 3 4], :);

    edge = [1 2; 2 3; 3 4; 4 1; 1 3; 2 4];
    for i = 1:6
        w(i) = integrate_tri(v1, edge(i,1), edge(i,2)) + integrate_tri(v2, edge(i,1), edge(i,2));
    end
    W = diag(w);
end

function [qp, w] = tri_gauss3(T)
    A = polyarea(T(:,1), T(:,2));

    bary = [1/6 1/6 2/3;
        1/6 2/3 1/6;
        2/3 1/6 1/6];

    qp = bary * T;
    w  = (A/3)*ones(3,1);
end

function val = integrate_tri(v, i, j)
    [qp, w] = tri_gauss3(v);

    val = 0;
    for k = 1:size(qp,1)
        [~, dlambdas] = wachspress2d(v, qp(k,:));

        dlambda1 = dlambdas(i,:);
        dlambda2 = dlambdas(j,:);

        val = val + w(k) * dot(dlambda1, dlambda2);
    end
end