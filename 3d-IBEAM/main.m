% Parameters
b        = [0; 0; 1];
k_values = 2:5;
a_values = (1/10).^(0:3);
nquad    = 2;

% Exact solution
u_true = @(x,y,z) sin(pi*x).*sin(pi*y).*sin(pi*z);
u_x    = @(x,y,z) pi*cos(pi*x).*sin(pi*y).*sin(pi*z);
u_y    = @(x,y,z) pi*sin(pi*x).*cos(pi*y).*sin(pi*z);
u_z    = @(x,y,z) pi*sin(pi*x).*sin(pi*y).*cos(pi*z);

% Error storage
H1_err   = zeros(length(k_values), length(a_values));
L2_err   = zeros(length(k_values), length(a_values));
h_values = zeros(length(k_values), 1);

%% Create output directory
paraview_dir = fullfile(pwd, 'Paraview');
if ~exist(paraview_dir, 'dir')
    mkdir(paraview_dir);
end

for idx = 1:length(k_values)
    k = k_values(idx);
    h = 2^(-k);
    h_values(idx) = h;

    [nodes, elements] = generate_mesh3d(h);
    n_nodes = size(nodes, 1);

    for a_idx = 1:length(a_values)
        a = a_values(a_idx);

        f_value = @(x) 3*a*pi^2*sin(pi*x(1))*sin(pi*x(2))*sin(pi*x(3)) ...
            - pi*(b(1)*cos(pi*x(1))*sin(pi*x(2))*sin(pi*x(3)) ...
                + b(2)*sin(pi*x(1))*cos(pi*x(2))*sin(pi*x(3)) ...
                + b(3)*sin(pi*x(1))*sin(pi*x(2))*cos(pi*x(3)));

        [u_h, ~] = ibeam3d(a, b, f_value, u_true, elements, nodes, n_nodes, nquad, false);

        H1_err(idx, a_idx) = H1_semi_error3d(u_h, u_x, u_y, u_z, elements, nodes, nquad);
        L2_err(idx, a_idx) = L2_error3d(u_h, u_true, elements, nodes, nquad);

        %% Export to ParaView
        u_ex_nodes = u_true(nodes(:,1), nodes(:,2), nodes(:,3));
        export_to_paraview(nodes, elements, u_h, u_ex_nodes, ...
            fullfile(paraview_dir, sprintf('ibeam3d_k%d_a%g.vtu', k, a)));
    end
end

% Convergence plots
plot_convergence(h_values, L2_err, 'L^2 Norm Error',      a_values, false);
plot_convergence(h_values, H1_err, 'H^1 Semi-norm Error', a_values, true);

%% Convergence plot
function plot_convergence(h_values, err, title_str, a_values, is_H1)
    figure('Name', title_str);
    loglog(h_values, err, '-o', 'LineWidth', 1.5); hold on
    if is_H1
        C = err(1,end) / h_values(1);
        loglog(h_values, C * h_values,    '--k', 'LineWidth', 2);
        ref_label = 'O(h)';
    else
        C = err(1,end) / h_values(1)^2;
        loglog(h_values, C * h_values.^2, '--r', 'LineWidth', 2);
        ref_label = 'O(h^2)';
    end
    xlabel('h'); ylabel('error'); title(title_str);
    legend_entries = arrayfun(@(a) sprintf('\\alpha = %g', a), a_values, 'UniformOutput', false);
    legend([legend_entries, {ref_label}], 'Location', 'best');
    grid on;
end