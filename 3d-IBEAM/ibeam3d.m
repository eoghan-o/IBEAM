function [u_ibeam, A_II] = ibeam3d(a, b, f_value, u_true, elements, nodes, n_nodes, nquad, use_poly)
    f = build_rhs(elements, nodes, n_nodes, f_value, nquad);
    A = assemble_stiffness(a, b, elements, nodes, n_nodes, nquad);
    if use_poly
        [interior, bdy_nodes, g] = classify_nodes_polymesher(nodes, n_nodes, u_true);
    else
        [interior, bdy_nodes, g] = classify_nodes(nodes, n_nodes, u_true);
    end
    A_II = A(interior, interior);
    A_IB = A(interior, bdy_nodes);
    f_I  = f(interior);
    u_I  = A_II \ (f_I - A_IB * g);
    u_ibeam = zeros(n_nodes, 1);
    u_ibeam(interior)  = u_I;
    u_ibeam(bdy_nodes) = g;
end

function [interior, bdy_nodes, g] = classify_nodes_polymesher(nodes, n_nodes, u_true)
    tol = 1e-3;
    bdy_nodes = find( ...
        nodes(:,1) < tol | nodes(:,1) > 1-tol | ...
        nodes(:,2) < tol | nodes(:,2) > 1-tol | ...
        nodes(:,3) < tol | nodes(:,3) > 1-tol);
    interior = setdiff(1:n_nodes, bdy_nodes);
    g = arrayfun(@(k) u_true(nodes(k,1), nodes(k,2), nodes(k,3)), bdy_nodes);
end

