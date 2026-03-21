function f = build_rhs(elements, nodes, n_nodes, f_value, n_quad)
    f = zeros(n_nodes,1);
    for T = 1:size(elements,1)
        indices = elements{T};
        v = nodes(indices, :);
    
        integrand = @(qp) wachspress2d(v, qp) * f_value(qp(1), qp(2));
        % Note: wachspress2d must return phi as a column vector here
    
        quadratures = polyintegrate(v, n_quad, integrand);
        f(indices) = f(indices) + quadratures;
    end
end

% function f = build_rhs_old(elements, nodes, n_nodes, f_value)
%     f = zeros(n_nodes,1);
% 
%     for T = 1:size(elements,1)
%         indices = elements(T, :);
%         vertices = nodes(indices, :);
% 
%         x_coords = vertices(:,1);
%         y_coords = vertices(:,2);
% 
%         x1 = min(x_coords);
%         x2 = max(x_coords);
%         y1 = min(y_coords);
%         y2 = max(y_coords);
% 
%         phi1 = @(x,y) f_value(x,y).*((x-x2).*(y-y2))./((x1-x2).*(y1-y2));
%         phi2 = @(x,y) f_value(x,y).*((x-x1).*(y-y2))./((x2-x1).*(y1-y2));
%         phi3 = @(x,y) f_value(x,y).*((x-x2).*(y-y1))./((x1-x2).*(y2-y1));
%         phi4 = @(x,y) f_value(x,y).*((x-x1).*(y-y1))./((x2-x1).*(y2-y1));
% 
%         f(indices(1)) = f(indices(1)) + integral2(phi1, x1, x2, y1, y2);
%         f(indices(2)) = f(indices(2)) + integral2(phi2, x1, x2, y1, y2);
%         f(indices(3)) = f(indices(3)) + integral2(phi3, x1, x2, y1 ,y2);
%         f(indices(4)) = f(indices(4)) + integral2(phi4, x1, x2, y1, y2);
%     end
% end