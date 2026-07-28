addpath('PolyMesher')
close all;
% set(groot, 'defaultFigureColor', 'white', ...
%            'defaultAxesColor', 'white', ...
%            'defaultTextColor', 'black', ...
%            'defaultAxesGridColor', [0.15 0.15 0.15], ...
%            'defaultAxesXColor', 'black', ...
%            'defaultAxesYColor', 'black')
% set(groot, 'defaultLegendColor', 'white', ...
%            'defaultLegendTextColor', 'black')

% Parameters
b          = [0; 1];
k_values   = 3:9;
a_values   = (1/10).^(0:3);

nquad      = 3;
use_poly   = true; %false;   % true: PolyMesher, false: structured quad mesh
perturb = 0.0;

% Exact solution
u_true = @(x,y) sin(pi*x).*sin(pi*y);
u_x    = @(x,y) pi*cos(pi*x).*sin(pi*y);
u_y    = @(x,y) pi*cos(pi*y).*sin(pi*x); 


% Error storage
H1_err = zeros(length(k_values), length(a_values));
L2_err = zeros(length(k_values), length(a_values));
h_values = zeros(length(k_values), 1);
for idx = 1:length(k_values)
    k = k_values(idx);
    
    if use_poly
        n_elements = 2^k;
        [elements, nodes, n_nodes, h_values(idx)] = generate_mesh(n_elements);
    else
        n = 2^k + 1;
        [elements, nodes, n_nodes] = generate_perturbed_mesh(n, perturb);
        h_values(idx) = 1 / (n - 1);
    end

    for a_idx = 1:length(a_values)
        a = a_values(a_idx);
    
        r2   = @(x,y) (x-0.5).^2 + (y-0.5).^2;
        k    = 2/sqrt(a);
        s    = @(x,y) k .* (0.0625 - r2(x,y));
        D    = @(x,y) 1 + s(x,y).^2;
        
        P    = @(x,y) 16.*x.*(1-x).*y.*(1-y);
        Px   = @(x,y) 16.*(1-2.*x).*y.*(1-y);
        Py   = @(x,y) 16.*x.*(1-x).*(1-2.*y);
        Pxx  = @(x,y) -32.*y.*(1-y);
        Pyy  = @(x,y) -32.*x.*(1-x);
        
        Q    = @(x,y) 0.5 + (1/pi).*atan(s(x,y));
        Qx   = @(x,y) (1/pi) .* (-2*k*(x-0.5)) ./ D(x,y);
        Qy   = @(x,y) (1/pi) .* (-2*k*(y-0.5)) ./ D(x,y);
        
        % d^2/dx^2 arctan(s) = [s_xx*(1+s^2) - 2s*s_x^2] / (1+s^2)^2
        % s_xx = s_yy = -2k  (constant!),  s_x = -2k(x-0.5),  s_y = -2k(y-0.5)
        Qxx  = @(x,y) (1/pi) .* ((-2*k).*D(x,y) - 2.*s(x,y).*(2*k*(x-0.5)).^2) ./ D(x,y).^2;
        Qyy  = @(x,y) (1/pi) .* ((-2*k).*D(x,y) - 2.*s(x,y).*(2*k*(y-0.5)).^2) ./ D(x,y).^2;
        
        %u_true = @(x,y) P(x,y) .* Q(x,y);
        %u_x    = @(x,y) Px(x,y).*Q(x,y)  + P(x,y).*Qx(x,y);
        %u_y    = @(x,y) Py(x,y).*Q(x,y)  + P(x,y).*Qy(x,y);
        
        uxx  = @(x,y) Pxx(x,y).*Q(x,y) + 2.*Px(x,y).*Qx(x,y) + P(x,y).*Qxx(x,y);
        uyy  = @(x,y) Pyy(x,y).*Q(x,y) + 2.*Py(x,y).*Qy(x,y) + P(x,y).*Qyy(x,y);
        
        %f_value = @(x,y) -a.*(uxx(x,y) + uyy(x,y)) + b(1).*u_x(x,y) + b(2).*u_y(x,y);
        
        f_value = @(x,y) 2*a*pi^2*sin(pi*x).*sin(pi*y) - pi*(b(1)*cos(pi*x).*sin(pi*y) + b(2)*sin(pi*x).*cos(pi*y));

        [u_ibeam, ~] = ibeam(a, b, f_value, u_true, elements, nodes, n_nodes, nquad, use_poly);
        H1_err(idx, a_idx) = H1_semi_error(u_ibeam, u_x, u_y, elements, nodes, nquad);
        L2_err(idx, a_idx) = L2_error(u_ibeam, u_true, elements, nodes, nquad);
    end
end

% Plot solution (structured mesh only)
if ~use_poly
    n = 2^k_values(end) + 1;
    [X, Y] = meshgrid(linspace(0, 1, n));
    figure('Name', 'Graph of Solution')
    surf(X, Y, reshape(u_ibeam, n, n));
    xlabel('x'); ylabel('y'); zlabel('u(x,y)');
    axis square; grid on
    title('Solution')

    % Plot (possibly perturbed) structured quad mesh
    % Elements are stored as triangle pairs — merge every two into a quad
    figure('Name', 'Mesh')
    hold on
    e = 1;
    while e < length(elements)
        t1 = elements{e};
        t2 = elements{e+1};
        % The quad corners are the union of both triangles' nodes
        corners = unique([t1, t2], 'stable');
        x_c = nodes(corners, 1);
        y_c = nodes(corners, 2);
        % Sort corners by angle around centroid to get correct polygon order
        cx = mean(x_c); cy = mean(y_c);
        angles = atan2(y_c - cy, x_c - cx);
        [~, order] = sort(angles);
        fill(x_c(order), y_c(order), 'w', 'EdgeColor', 'k', 'LineWidth', 0.5);
        e = e + 2;
    end
    xlabel('x'); ylabel('y');
    axis equal tight; grid off
    title(sprintf('Mesh (n = %d, perturb = %g)', n, perturb))
    hold off
else
    % Unstructured mesh — plot solution
    figure('Name', 'Graph of Solution')
    patch('Faces', cell2mat(cellfun(@(e) [e, nan(1, max(cellfun(@length, elements)) - length(e))], elements, 'UniformOutput', false)), ...
        'Vertices', nodes, ...
        'FaceVertexCData', u_ibeam, ...
        'FaceColor', 'interp', ...
        'EdgeColor', 'k');
    colorbar
    xlabel('x'); ylabel('y');
    axis square; grid on
    title('Solution')

    % Plot unstructured poly mesh
    figure('Name', 'Mesh')
    max_verts = max(cellfun(@length, elements));
    face_matrix = cell2mat(cellfun(@(e) [e, nan(1, max_verts - length(e))], elements, 'UniformOutput', false));
    patch('Faces', face_matrix, ...
        'Vertices', nodes, ...
        'FaceColor', 'white', ...
        'EdgeColor', 'k', ...
        'LineWidth', 0.5);
    xlabel('x'); ylabel('y');
    axis square; grid on
    title(sprintf('Polygonal Mesh (%d elements)', length(elements)))
end

% Convergence plots
if use_poly
    % Build vector of element counts across k_values
    n_el_vec = arrayfun(@(k) 2^k, k_values)';
    plot_convergence_nel(n_el_vec, L2_err, 'L^2 Norm Error', a_values, false);
    plot_convergence_nel(n_el_vec, H1_err, 'H^1 Semi-norm Error', a_values, true);
else
    plot_convergence(h_values, L2_err, 'L^2 Norm Error', a_values, false);
    plot_convergence(h_values, H1_err, 'H^1 Semi-norm Error', a_values, true);
end

function plot_convergence(h_values, err, title_str, a_values, is_H1)
    figure('Name', title_str);
    loglog(h_values, err, '-o', 'LineWidth', 1.5); hold on
    if is_H1
        C = err(end, end) / h_values(end);          % anchor at finest h
        loglog(h_values, C * h_values,    '--k', 'LineWidth', 2);
        ref_label = 'O(h)';
    else
        C = err(end, end) / h_values(end)^2;        % anchor at finest h
        loglog(h_values, C * h_values.^2, '--r', 'LineWidth', 2);
        ref_label = 'O(h^2)';
    end
    xlabel('h'); ylabel('error'); title(title_str);
    legend_entries = arrayfun(@(a) sprintf('\\alpha = %g', a), a_values, 'UniformOutput', false);
    legend([legend_entries, {ref_label}], 'Location', 'best');
    grid on;
end

function plot_convergence_nel(n_el_vec, err, title_str, a_values, is_H1)
    figure('Name', title_str);
    loglog(n_el_vec, err, '-o', 'LineWidth', 1.5); hold on
    if is_H1
        C = err(end, end) / n_el_vec(end)^(-0.5);
        loglog(n_el_vec, C * n_el_vec.^(-0.5), '--k', 'LineWidth', 2);
        ref_label = 'O(N^{-1/2})';
    else
        C = err(end, end) / n_el_vec(end)^(-1);
        loglog(n_el_vec, C * n_el_vec.^(-1), '--r', 'LineWidth', 2);
        ref_label = 'O(N^{-1})';
    end
    xlabel('Number of elements N'); ylabel('error'); title(title_str);
    legend_entries = arrayfun(@(a) sprintf('\\alpha = %g', a), a_values, 'UniformOutput', false);
    legend([legend_entries, {ref_label}], 'Location', 'best');
    grid on;
end

%save('squares.mat',   'h_values', 'L2_err', 'H1_err', 'a_values');  % use_poly=false, perturb=0
%save('quads.mat',     'h_values', 'L2_err', 'H1_err', 'a_values');  % use_poly=false, perturb>0
save('polygons.mat',  'h_values', 'L2_err', 'H1_err', 'a_values');  % use_poly=true
