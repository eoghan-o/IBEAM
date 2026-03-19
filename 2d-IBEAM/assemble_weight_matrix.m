function W = assemble_weight_matrix(v, n)
    w = zeros(6,1);

    v1 = v([1 2 3], :);
    v2 = v([1 3 4], :);

    [qp1, w1] = triquad(v1, n);
    [qp2, w2] = triquad(v2, n);

    edge = [1 2; 2 3; 3 4; 4 1; 1 3; 2 4];
    for k = 1:6
        i = edge(k,1);
        j = edge(k,2);
        
        I1 = 0;
        I2 = 0;
        for p = 1:n^2
            [~, dphi1] = wachspress2d(v, qp1(p,:));
            [~, dphi2] = wachspress2d(v, qp2(p,:));
    
            I1 = I1 + w1(p)*dot(dphi1(i,:), dphi1(j,:));
            I2 = I2 + w2(p)*dot(dphi2(i,:), dphi2(j,:));
        end
        w(k) = -1*(I1 + I2);
    end
    W = diag(w);
end