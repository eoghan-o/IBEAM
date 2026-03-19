function [outputArg1,outputArg2] = L2_error(u, u_true, elements, nodes, n_quad)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    total = 0;
    for T = 1:size(elements,1)
        indices = elements(T, :);
        v = nodes(indices, :);
        u_local = u(indices);
        
        v1 = v([1 2 3], :);
        v2 = v([1 3 4], :);

        [qp1, w1] = triquad(v1, n_quad);
        [qp2, w2] = triquad(v2, n_quad);

        for p = 1:n_quad^2
            [phi1, ~] = wachspress2d(v, qp1(p,:));
            [phi2, ~] = wachspress2d(v, qp2(p,:));

            u1_local = u1 .* phi1;
        end

    end
end