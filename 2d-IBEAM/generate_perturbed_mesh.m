function [elements, nodes, n_nodes] = generate_perturbed_mesh(n, perturb)
    [X,Y] = meshgrid(linspace(0,1,n));
    X = X'; Y = Y';
    nodes = [X(:), Y(:)];

    for j = 2:n-1
        for i = 2:n-1
            idx = i + (j-1)*n;
            nodes(idx, 1) = nodes(idx, 1) + (rand - 0.5)*perturb;
            nodes(idx, 2) = nodes(idx, 2) + (rand - 0.5)*perturb;
        end
    end
    n_nodes = size(nodes,1);
    
    n_elements = (n-1)^2;
    elements = zeros(n_elements,4);
    e = 1;
    for j = 1:n-1
        for i = 1:n-1
            n1 = i + (j-1)*n;
            n2 = i + 1 + (j-1)*n;
            n3 = i + 1 + j*n;
            n4 = i + j*n;
            elements(e,:) = [n1 n2 n3 n4];
            e = e + 1;
        end
    end
end