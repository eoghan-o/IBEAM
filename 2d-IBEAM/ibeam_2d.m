addpath('PolyMesher')
close all;

% Parameters
b          = [0; 1];
k_values   = 3:9;
a_values   = (1/10).^(0:2);
nquad      = 3;
use_poly   = true; %false;   % true: PolyMesher, false: structured quad mesh

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
        [elements, nodes, n_nodes] = generate_perturbed_mesh(n, 0);
        h_values(idx) = 1 / (n - 1);
    end

    for a_idx = 1:length(a_values)
        a = a_values(a_idx);

        f_value = @(x,y) 2*a*pi^2*sin(pi*x).*sin(pi*y) - ...
            pi*(b(1)*cos(pi*x).*sin(pi*y) + b(2)*sin(pi*x).*cos(pi*y));

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
else
    % Unstructured mesh
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
end

% Convergence plots
plot_convergence(h_values, L2_err, 'L^2 Norm Error', a_values);
plot_convergence(h_values, H1_err, 'H^1 Semi-norm Error', a_values);

function plot_convergence(h_values, err, title_str, a_values)
    C1 = err(1,1) / h_values(1);
    C2 = err(1,1) / h_values(1)^2;
    figure('Name', title_str)
    loglog(h_values, err, 'LineWidth', 1.5); hold on
    loglog(h_values, C1*h_values,    '--k', 'LineWidth', 1.2)
    loglog(h_values, C2*h_values.^2, '--r', 'LineWidth', 1.2)
    xlabel('h'); ylabel('error'); title(title_str)
    legend_entries = arrayfun(@(a) sprintf('\\alpha = %g', a), a_values, 'UniformOutput', false);
    legend([legend_entries, {'O(h)', 'O(h^2)'}], 'Location', 'best')
    grid on
end