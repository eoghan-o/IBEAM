close all;

b = [1; 2];
f_value = @(x,y) 1 + (x-x);
k_values = 4:6;
a_values = (1/10).^(2:5);

num_plots = length(k_values);
H1_ibeam_k = zeros(length(k_values),length(a_values));
L2_ibeam_k = zeros(length(k_values),length(a_values));

for idx = 1:num_plots
    k = k_values(idx);
    n = 2^k + 1;
    h = 1/(n-1);
    perturb = .5*h;
    [elements, nodes, n_nodes] = generate_perturbed_mesh(n, perturb);

    for a_idx = 1:length(a_values)
        a = a_values(a_idx);
        f_value = @(x,y) 0 + a*2*sin(x+y) - (b(1) + b(2)) * cos(x + y);
        [u_ibeam, A] = ibeam(k, a, b, f_value, elements, nodes, n_nodes);
        H1_ibeam_k(idx,a_idx) = H1_error(u_ibeam,elements,nodes);
        L2_ibeam_k(idx,a_idx) = L2_error(u_ibeam,elements,nodes);
    end
end
[X, Y] = meshgrid(linspace(0, 1, n));
u = ones(n*n,1);

figure('name', 'Graph of Solution')
surf(X, Y, reshape(u_ibeam, n, n));  hold on
xlabel('x'); ylabel('y'); zlabel('u(x,y)'); 
axis square
grid on

h_values = 1 ./ (2.^k_values);
C1 = L2_ibeam_k(1,1) / h_values(1);
C2 = L2_ibeam_k(1,1) / h_values(1)^2;

figure('Name', 'L^2 norm error')
loglog(h_values, L2_ibeam_k(:,1), 'LineWidth', 1.5); hold on
loglog(h_values, L2_ibeam_k(:,2), 'LineWidth', 1.5)
loglog(h_values, L2_ibeam_k(:,3), 'LineWidth', 1.5)
loglog(h_values, L2_ibeam_k(:,4), 'LineWidth', 1.5)
loglog(h_values, C1*h_values, '--k', 'LineWidth', 1.2)
loglog(h_values, C2*h_values.^2, '--r', 'LineWidth', 1.2)
xlabel('h'); ylabel('error'); title('L^2 norm error')
legend('\alpha = 0.1', '\alpha = 0.01', '\alpha = 1E-3', '\alpha = 1E-4', 'O(h)', 'O(h^2)', 'Location', 'best'); grid on

function value = true_solution(a, x, y)
    %value = x.*(1 - exp((y-1)./a))./(1-exp(-2/a));
    %value = 0;
    value = sin(x+y);
end

function err = L2_error(u,elements,nodes)
    total = 0;
    for T = 1:size(elements,1)
        indices = elements(T, :);
        vertices = nodes(indices, :);
        u_local = u(indices);

        x_coords = vertices(:,1);
        y_coords = vertices(:,2);

        x1 = min(x_coords);
        x2 = max(x_coords);
        y1 = min(y_coords);
        y2 = max(y_coords);

        phi1 = @(x,y) ((x-x2).*(y-y2))./((x1-x2).*(y1-y2));
        phi2 = @(x,y) ((x-x1).*(y-y2))./((x2-x1).*(y1-y2));
        phi3 = @(x,y) ((x-x2).*(y-y1))./((x1-x2).*(y2-y1));
        phi4 = @(x,y) ((x-x1).*(y-y1))./((x2-x1).*(y2-y1));
        u_true = @(x,y) sin(x+y);

        uh = @(x,y) u_local(1)*phi1(x,y) + ...
            u_local(2)*phi2(x,y) + ...
            u_local(3)*phi3(x,y) + ...
            u_local(4)*phi4(x,y);

        integrand = @(x,y) (uh(x,y) - u_true(x,y)).^2;

        total = total + integral2(integrand, x1, x2, y1, y2);
    end
    err = sqrt(total);
end

function err = H1_error(u,elements,nodes)
    total = 0;
    for T = 1:size(elements,1)
        indices = elements(T, :);
        vertices = nodes(indices, :);
        u_local = u(indices);

        x_coords = vertices(:,1);
        y_coords = vertices(:,2);

        x1 = min(x_coords);
        x2 = max(x_coords);
        y1 = min(y_coords);
        y2 = max(y_coords);

        phi1 = @(x,y) ((x-x2).*(y-y2))./((x1-x2).*(y1-y2));
        phi2 = @(x,y) ((x-x1).*(y-y2))./((x2-x1).*(y1-y2));
        phi3 = @(x,y) ((x-x2).*(y-y1))./((x1-x2).*(y2-y1));
        phi4 = @(x,y) ((x-x1).*(y-y1))./((x2-x1).*(y2-y1));

        phi1_x = @(x,y) (y-y2)./((x1-x2).*(y1-y2));
        phi1_y = @(x,y) (x-x2)./((x1-x2).*(y1-y2));
        phi2_x = @(x,y) (y-y2)./((x2-x1).*(y1-y2));
        phi2_y = @(x,y) (x-x1)./((x2-x1).*(y1-y2));
        phi3_x = @(x,y) (y-y1)./((x1-x2).*(y2-y1));
        phi3_y = @(x,y) (x-x2)./((x1-x2).*(y2-y1));
        phi4_x = @(x,y) (y-y1)./((x2-x1).*(y2-y1));
        phi4_y = @(x,y) (x-x1)./((x2-x1).*(y2-y1));

        uh = @(x,y) u_local(1)*phi1(x,y) + ...
            u_local(2)*phi2(x,y) + ...
            u_local(3)*phi3(x,y) + ...
            u_local(4)*phi4(x,y);

        grad_uh_x = @(x,y) u_local(1)*phi1_x(x,y) + ...
            u_local(2)*phi2_x(x,y) + ...
            u_local(3)*phi3_x(x,y) + ...
            u_local(4)*phi4_x(x,y);

        grad_uh_y = @(x,y) u_local(1)*phi1_y(x,y) + ...
            u_local(2)*phi2_y(x,y) + ...
            u_local(3)*phi3_y(x,y) + ...
            u_local(4)*phi4_y(x,y);

        u_true = @(x,y) sin(x + y);
        grad_u_x = @(x,y) cos(x+y);
        grad_u_y = @(x,y) cos(x+y);

        % u_true = @(x,y) (x-x) + 0;
        % grad_u_x = @(x,y) (x-x) + 0;
        % grad_u_y = @(x,y) (x-x) + 0;

        integrand = @(x,y) (uh(x,y) - u_true(x,y)).^2 + ...
            (grad_uh_x(x,y) - grad_u_x(x,y)).^2 + ...
            (grad_uh_y(x,y) - grad_u_y(x,y)).^2;

        total = total + integral2(integrand, x1, x2, y1, y2);
    end
    err = sqrt(total);
end

function [elements, nodes, n_nodes] = generate_mesh(n)
    [X,Y] = meshgrid(linspace(0,1,n));
    nodes = [Y(:), X(:)];
    n_nodes = size(nodes,1);

    elements = [];
    for j = 1:n-1
        for i = 1:n-1
            n1 = (j-1)*n + i;
            n2 = n1 + 1;
            n3 = n2 + n;
            n4 = n1 + n;

            elements = [elements; n1 n2 n3 n4];
        end
    end
end

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

function B = bernoulli(p)
    if abs(p) < 10^(-8)
        B = 1 - p/2 + p^2;
    else
        B = p/(exp(p) - 1);
    end
end

function f = build_rhs(n, f_value)
    [elements, nodes, n_nodes] = generate_mesh(n);

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

function [phi, dphi] = wachspress2d(v,x)
% Written by M. Floater, et al.
% Evaluate Wachspress basis functions and their gradients in a convex
% polygon
%
% Inputs:
% v     : [x1 y1; x2 y2; ...; xn yn], the n vertices of the polygon in ccw
% x     : [x(1) x(2)], the point at which the basis functions are computed
% Outputs:
% phi   : output basis functions = [phi_1; ...; phi_n]
% dphi  : output gradient of basis functions = [dphi_1; ...; dphi_n]
    n = size(v,1);
    w = zeros(n,1);
    R = zeros(n,2);
    phi = zeros(n,1);
    dphi = zeros(n,2);

    un = getNormals(v);
    p = zeros(n,2);
    for i = 1:n
        h = dot(v(i,:) - x,un(i,:));
        p(i,:) = un(i,:) / h;
    end

    for i = 1:n
        im1 = mod(i-2,n) + 1;
        w(i) = det([p(im1,:); p(i,:)]);
        R(i,:) = p(im1,:) + p(i,:);
    end

    wsum = sum(w);
    phi = w/wsum;

    phiR = phi' * R;
    for k = 1:2
        dphi(:,k) = phi .* (R(:,k) - phiR(:,k));
    end
end

function un = getNormals(v)
% Written by M. Floater, et al.
% Function to compute the outword unit normal to each edge
    n = size(v,1);
    un = zeros(n,2);
    for i = 1:n
        d = v(mod(i,n) + 1,:) - v(i,:);
        un(i,:) = [d(2) -d(1)]/norm(d);
    end
end

function [qp, w] = tri_gauss3(T)
    A = polyarea(T(:,1), T(:,2));

    bary = [1/6 1/6 2/3;
        1/6 2/3 1/6;
        2/3 1/6 1/6];

    qp = bary * T;
    w  = (A/3)*ones(3,1);
end

function val = integrate_tri(v, i, j)
    [qp, w] = tri_gauss3(v);

    val = 0;
    for k = 1:size(qp,1)
        [~, dlambdas] = wachspress2d(v, qp(k,:));

        dlambda1 = dlambdas(i,:);
        dlambda2 = dlambdas(j,:);

        val = val + w(k) * dot(dlambda1, dlambda2);
    end
end

function W = compute_w(vertices)
    w = zeros(6,1);

    v1 = vertices([1 2 3], :);
    v2 = vertices([1 3 4], :);

    edge = [1 2; 2 3; 3 4; 4 1; 1 3; 2 4];
    for i = 1:6
        w(i) = integrate_tri(v1, edge(i,1), edge(i,2)) + integrate_tri(v2, edge(i,1), edge(i,2));
    end
    W = diag(w);
end

function [u_ibeam, A] = ibeam(k, a, b, f_value, elements, nodes, n_nodes)
    n = 2^k + 1;
    %f = zeros(n_nodes,1);
    f = build_rhs(n, f_value);

    A = sparse(n_nodes, n_nodes);

    %G_T = [-1 1 0 0; -1 0 1 0; 0 -1 0 1; 0 0 -1 1; -1 0 0 1; 0 -1 1 0];
    G_T = [-1 1 0 0; 0 -1 1 0; 0 0 -1 1; 1 0 0 -1; -1 0 1 0; 0 -1 0 1];
    local_edge = [2 1; 3 2; 4 3; 1 4; 3 1; 2 4];
    w = [1/6 1/6 1/6 1/6 1/3 1/3]';
    W = diag(w);

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

        % w = compute_W(vertices, quad_pts, quad_w);
        A_T = a*G_T'*W*Gb_T;
        A(indices, indices) = A(indices, indices) + A_T;
    end

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

    rhs_reduced = f_I - A_IB * g;
    u_I = A_II \ rhs_reduced;

    u_ibeam = zeros(n_nodes,1);
    u_ibeam(interior) = u_I;
    u_ibeam(bdy_nodes) = g;
end