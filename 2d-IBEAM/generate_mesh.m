function [elements, nodes, n_nodes, h] = generate_mesh(n_elements)
    [~, nodes, elements] = evalc('PolyMesher(@square_domain, n_elements, 100)');
    
    % Compute average edge length as h
    total_edge = 0;
    n_edges = 0;
    for T = 1:length(elements)
        idx = elements{T};
        n_v = length(idx);
        for i = 1:n_v
            j = mod(i, n_v) + 1;
            total_edge = total_edge + norm(nodes(idx(j),:) - nodes(idx(i),:));
            n_edges = n_edges + 1;
        end
    end
    h = total_edge / n_edges;
    
    % Enforce CCW orientation
    for T = 1:length(elements)
        idx = elements{T};
        v = nodes(idx,:);
        n_v = size(v,1);
        area = 0;
        for i = 1:n_v
            j = mod(i, n_v) + 1;
            area = area + v(i,1)*v(j,2) - v(j,1)*v(i,2);
        end
        if area < 0
            elements{T} = fliplr(idx);
        end
    end
    
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