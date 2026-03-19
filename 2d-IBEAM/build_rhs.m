function f = build_rhs(elements, nodes, n_nodes, f_value, n_quad)
    f = zeros(n_nodes,1);

    for T = 1:size(elements,1)
        indices = elements(T, :);
        v = nodes(indices, :);

        v1 = v([1 2 3], :);
        v2 = v([1 3 4], :);

        [qp1, w1] = triquad(v1, n_quad);
        [qp2, w2] = triquad(v2, n_quad);
        
        quadratures = zeros(size(elements(T,:)',1),1);
        for p = 1:n_quad^2
            [phi1, ~] = wachspress2d(v, qp1(p,:));
            [phi2, ~] = wachspress2d(v, qp2(p,:));
            f1 = f_value(qp1(p,1), qp1(p,2));
            f2 = f_value(qp2(p,1), qp2(p,2));

            quadratures = quadratures + w1(p)*f1*phi1 + w2(p)*f2*phi2;
        end
        f(indices) = f(indices) + quadratures;
    end
    %f_old = build_rhs_old(elements, nodes, n_nodes, f_value);
end

function f = build_rhs_old(elements, nodes, n_nodes, f_value)
    f = zeros(n_nodes,1);

    for T = 1:size(elements,1)
        indices = elements(T, :);
        vertices = nodes(indices, :);

        x_coords = vertices(:,1);
        y_coords = vertices(:,2);

        x1 = min(x_coords);
        x2 = max(x_coords);
        y1 = min(y_coords);
        y2 = max(y_coords);

        phi1 = @(x,y) f_value(x,y).*((x-x2).*(y-y2))./((x1-x2).*(y1-y2));
        phi2 = @(x,y) f_value(x,y).*((x-x1).*(y-y2))./((x2-x1).*(y1-y2));
        phi3 = @(x,y) f_value(x,y).*((x-x2).*(y-y1))./((x1-x2).*(y2-y1));
        phi4 = @(x,y) f_value(x,y).*((x-x1).*(y-y1))./((x2-x1).*(y2-y1));

        f(indices(1)) = f(indices(1)) + integral2(phi1, x1, x2, y1, y2);
        f(indices(2)) = f(indices(2)) + integral2(phi2, x1, x2, y1, y2);
        f(indices(3)) = f(indices(3)) + integral2(phi3, x1, x2, y1 ,y2);
        f(indices(4)) = f(indices(4)) + integral2(phi4, x1, x2, y1, y2);
    end
end