close all;

% indices = elements(1, :);
% v = nodes(indices, :);
% 
% i = 1;
% j = 2;
% 
% t1 = v([1 2 3], :);
% t2 = v([1 3 4], :);
% result = zeros(6,1);
% 
% for n = 2:7
%     [qp1, w1] = triquad(n, t1);
%     [qp2, w2] = triquad(n, t2);
% 
%     I1 = 0; I2 = 0;
%     for k = 1:n^2
%         [phi1, dphi1] = wachspress2d(v, qp1(k,:));
%         [phi2, dphi2] = wachspress2d(v, qp2(k,:));
% 
%         I1 = I1 + w1(k)*dot(dphi1(i,:), dphi1(j,:));
%         I2 = I2 + w2(k)*dot(dphi2(i,:), dphi2(j,:));
%     end
%     result(n-1) = I1 + I2;
% end

b = [1; 2];
k_values = 4:6;
a_values = (1/10).^(2:3);

nquad = 3;

num_plots = length(k_values);
H1_ibeam_k = zeros(length(k_values),length(a_values));
L2_ibeam_k = zeros(length(k_values),length(a_values));

for idx = 1:num_plots
    k = k_values(idx);
    n = 2^k + 1;
    h = 1/(n-1);
    perturb = 0;
    [elements, nodes, n_nodes] = generate_perturbed_mesh(n, perturb);

    for a_idx = 1:length(a_values)
        a = a_values(a_idx);
        %f_value = @(x,y) (x-x) + 1; 
        f_value = @(x,y) 0 + a*2*sin(x+y) - (b(1) + b(2)) * cos(x + y);
        [u_ibeam, A] = ibeam(k, a, b, f_value, elements, nodes, n_nodes, nquad);

        [X,Y] = meshgrid(linspace(0,1,n));
        figure('Name', sprintf('Graph of Solution: n=%d, alpha=%.3f', n, a))
        surf(X,Y, reshape(u_ibeam, n, n)); hold on
        xlabel('x'); ylabel('y'); zlabel('u(x,y)');
        axis square
        grid on
        title(sprintf('Solution for n = %d, \\alpha = %.3f', n, a))
        
        % H1_ibeam_k(idx,a_idx) = H1_error(u_ibeam,elements,nodes);
        % L2_ibeam_k(idx,a_idx) = L2_error(u_ibeam,elements,nodes);
    end
end
[X, Y] = meshgrid(linspace(0, 1, n));
u = ones(n*n,1);

figure('name', 'Graph of Solution')
surf(X, Y, reshape(u_ibeam, n, n));  hold on
xlabel('x'); ylabel('y'); zlabel('u(x,y)'); 
axis square
grid on

% h_values = 1 ./ (2.^k_values);
% C1 = L2_ibeam_k(1,1) / h_values(1);
% C2 = L2_ibeam_k(1,1) / h_values(1)^2;

% figure('Name', 'L^2 norm error')
% loglog(h_values, L2_ibeam_k(:,1), 'LineWidth', 1.5); hold on
% loglog(h_values, L2_ibeam_k(:,2), 'LineWidth', 1.5)
% loglog(h_values, L2_ibeam_k(:,3), 'LineWidth', 1.5)
% loglog(h_values, L2_ibeam_k(:,4), 'LineWidth', 1.5)
% loglog(h_values, C1*h_values, '--k', 'LineWidth', 1.2)
% loglog(h_values, C2*h_values.^2, '--r', 'LineWidth', 1.2)
% xlabel('h'); ylabel('error'); title('L^2 norm error')
% legend('\alpha = 0.1', '\alpha = 0.01', '\alpha = 1E-3', '\alpha = 1E-4', 'O(h)', 'O(h^2)', 'Location', 'best'); grid on

function value = true_solution(a, x, y)
    %value = x.*(1 - exp((y-1)./a))./(1-exp(-2/a));
    %value = 0;
    %value = 1;
    value = sin(x + y);
end

% function err = L2_error(u,elements,nodes)
%     total = 0;
%     for T = 1:size(elements,1)
%         indices = elements(T, :);
%         vertices = nodes(indices, :);
%         u_local = u(indices);
% 
%         x_coords = vertices(:,1);
%         y_coords = vertices(:,2);
% 
%         x1 = min(x_coords);
%         x2 = max(x_coords);
%         y1 = min(y_coords);
%         y2 = max(y_coords);
% 
%         phi1 = @(x,y) ((x-x2).*(y-y2))./((x1-x2).*(y1-y2));
%         phi2 = @(x,y) ((x-x1).*(y-y2))./((x2-x1).*(y1-y2));
%         phi3 = @(x,y) ((x-x2).*(y-y1))./((x1-x2).*(y2-y1));
%         phi4 = @(x,y) ((x-x1).*(y-y1))./((x2-x1).*(y2-y1));
%         u_true = @(x,y) sin(x+y);
% 
%         uh = @(x,y) u_local(1)*phi1(x,y) + ...
%             u_local(2)*phi2(x,y) + ...
%             u_local(3)*phi3(x,y) + ...
%             u_local(4)*phi4(x,y);
% 
%         integrand = @(x,y) (uh(x,y) - u_true(x,y)).^2;
% 
%         total = total + integral2(integrand, x1, x2, y1, y2);
%     end
%     err = sqrt(total);
% end
% 
% function err = H1_error(u,elements,nodes)
%     total = 0;
%     for T = 1:size(elements,1)
%         indices = elements(T, :);
%         vertices = nodes(indices, :);
%         u_local = u(indices);
% 
%         x_coords = vertices(:,1);
%         y_coords = vertices(:,2);
% 
%         x1 = min(x_coords);
%         x2 = max(x_coords);
%         y1 = min(y_coords);
%         y2 = max(y_coords);
% 
%         phi1 = @(x,y) ((x-x2).*(y-y2))./((x1-x2).*(y1-y2));
%         phi2 = @(x,y) ((x-x1).*(y-y2))./((x2-x1).*(y1-y2));
%         phi3 = @(x,y) ((x-x2).*(y-y1))./((x1-x2).*(y2-y1));
%         phi4 = @(x,y) ((x-x1).*(y-y1))./((x2-x1).*(y2-y1));
% 
%         phi1_x = @(x,y) (y-y2)./((x1-x2).*(y1-y2));
%         phi1_y = @(x,y) (x-x2)./((x1-x2).*(y1-y2));
%         phi2_x = @(x,y) (y-y2)./((x2-x1).*(y1-y2));
%         phi2_y = @(x,y) (x-x1)./((x2-x1).*(y1-y2));
%         phi3_x = @(x,y) (y-y1)./((x1-x2).*(y2-y1));
%         phi3_y = @(x,y) (x-x2)./((x1-x2).*(y2-y1));
%         phi4_x = @(x,y) (y-y1)./((x2-x1).*(y2-y1));
%         phi4_y = @(x,y) (x-x1)./((x2-x1).*(y2-y1));
% 
%         uh = @(x,y) u_local(1)*phi1(x,y) + ...
%             u_local(2)*phi2(x,y) + ...
%             u_local(3)*phi3(x,y) + ...
%             u_local(4)*phi4(x,y);
% 
%         grad_uh_x = @(x,y) u_local(1)*phi1_x(x,y) + ...
%             u_local(2)*phi2_x(x,y) + ...
%             u_local(3)*phi3_x(x,y) + ...
%             u_local(4)*phi4_x(x,y);
% 
%         grad_uh_y = @(x,y) u_local(1)*phi1_y(x,y) + ...
%             u_local(2)*phi2_y(x,y) + ...
%             u_local(3)*phi3_y(x,y) + ...
%             u_local(4)*phi4_y(x,y);
% 
%         u_true = @(x,y) sin(x + y);
%         grad_u_x = @(x,y) cos(x+y);
%         grad_u_y = @(x,y) cos(x+y);
% 
%         % u_true = @(x,y) (x-x) + 0;
%         % grad_u_x = @(x,y) (x-x) + 0;
%         % grad_u_y = @(x,y) (x-x) + 0;
% 
%         integrand = @(x,y) (uh(x,y) - u_true(x,y)).^2 + ...
%             (grad_uh_x(x,y) - grad_u_x(x,y)).^2 + ...
%             (grad_uh_y(x,y) - grad_u_y(x,y)).^2;
% 
%         total = total + integral2(integrand, x1, x2, y1, y2);
%     end
%     err = sqrt(total);
% end

function B = bernoulli(p)
    if abs(p) < 10^(-8)
        B = 1 - p/2 + p^2;
    else
        B = p/(exp(p) - 1);
    end
end

function [u_ibeam, A] = ibeam(k, a, b, f_value, elements, nodes, n_nodes, nquad)
    n = 2^k + 1;
    %f = zeros(n_nodes,1);
    [f, f_old] = build_rhs(elements, nodes, n_nodes, f_value, nquad);
    norm(f - f_old)
    
    A = sparse(n_nodes, n_nodes);
    A_old = sparse(n_nodes, n_nodes);

    %G_T = [-1 1 0 0; -1 0 1 0; 0 -1 0 1; 0 0 -1 1; -1 0 0 1; 0 -1 1 0];
    G_T = [-1 1 0 0; 0 -1 1 0; 0 0 -1 1; 1 0 0 -1; -1 0 1 0; 0 -1 0 1];
    local_edge = [2 1; 3 2; 4 3; 1 4; 3 1; 2 4];
    w_old = [1/6 1/6 1/6 1/6 1/3 1/3]';
    W_old = diag(w_old);

    for T = 1:size(elements, 1)
        indices = elements(T, :);
        vertices = nodes(indices, :);
        Gb_T = zeros(6,4);
        for e = 1:6
            e_vec = (vertices(local_edge(e,1),:)-vertices(local_edge(e,2),:))';
            p_e = a^(-1)* b'*(e_vec);
            Gb_T(e,local_edge(e,1)) = bernoulli(-p_e);
            Gb_T(e,local_edge(e,2)) = -bernoulli(p_e);
        end
        
        W = assemble_weight_matrix(vertices, nquad);
        A_T = a*G_T'*W*Gb_T;
        A_T_old = a*G_T'*W_old*Gb_T;
        A(indices, indices) = A(indices, indices) + A_T;
        A_old(indices, indices) = A_old(indices, indices) + A_T_old;
    end
    % norm(f - f_old)
    % A_diff = full(A - A_old);
    % norm(A_diff)

    tol = 1e-12;
    bdy_nodes = find(...
        abs(nodes(:,1) - 0) < tol | abs(nodes(:,1) - 1) < tol | ...
        abs(nodes(:,2) - 0) < tol | abs(nodes(:,2) - 1) < tol);
    interior = setdiff(1:n_nodes, bdy_nodes);

    g = zeros(size(bdy_nodes));
    for k = 1:length(bdy_nodes)
        idx = bdy_nodes(k);
        x = nodes(idx, 1);
        y = nodes(idx, 2);
        g(k) = true_solution(a, x, y);
    end

    A_II = A(interior, interior);
    A_IB = A(interior, bdy_nodes);
    f_I = f(interior);

    A_II_old = A_old(interior, interior);
    A_IB_old = A(interior, bdy_nodes);
    f_I_old = f_old(interior);

    rhs_reduced = f_I - A_IB * g;
    u_I = A_II \ rhs_reduced;

    rhs_reduced_old = f_I_old - A_IB_old * g;
    u_I_old = A_II_old \ rhs_reduced_old;

    u_ibeam = zeros(n_nodes,1);
    u_ibeam(interior) = u_I;
    u_ibeam(bdy_nodes) = g;

    u_ibeam_old = zeros(n_nodes,1);
    u_ibeam_old(interior) = u_I_old;
    u_ibeam_old(bdy_nodes) = g;

    norm(u_ibeam - u_ibeam_old)
end