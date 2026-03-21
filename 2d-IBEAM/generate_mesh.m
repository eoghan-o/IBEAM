function [elements, nodes, n_nodes] = generate_mesh(n)
    n_elements = (n-1)^2;
    [nodes, elements] = PolyMesher(@SquareDomain, n_elements, 100);
    n_nodes = size(nodes, 1);
end