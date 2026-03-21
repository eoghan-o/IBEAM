function [elements, nodes, n_nodes] = generate_mesh(n_elements)
    [captured_output, nodes, elements] = evalc('PolyMesher(@square_domain, n_elements, 2000)');
    tol = 1e-10;
    nodes(abs(nodes(:,1)) < tol, 1) = 0;
    nodes(abs(nodes(:,1) - 1) < tol, 1) = 1;
    nodes(abs(nodes(:,2)) < tol, 2) = 0;
    nodes(abs(nodes(:,2) - 1) < tol, 2) = 1;
    n_nodes = size(nodes, 1);
end

function x = square_domain(Demand, Arg)
    BdBox = [0 1 0 1];
    switch Demand
        case 'Dist';  x = dRectangle(Arg, BdBox(1), BdBox(2), BdBox(3), BdBox(4));
        case 'BC';    x = {[], []};
        case 'BdBox'; x = BdBox;
        case 'PFix';  x = [];
    end
end
