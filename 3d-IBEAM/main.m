% Parameters
b        = [0; 0; 1];
k_values = 2:3;
a_values = (1/10).^(1:1);
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
    end
end

% Convergence plots
plot_convergence(h_values, L2_err, 'L^2 Norm Error',       a_values, false);
plot_convergence(h_values, H1_err, 'H^1 Semi-norm Error',  a_values, true);

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

%% Plot solution on final mesh
k   = k_values(end);
h   = 1/k;
a   = a_values(end);

[nodes, elements] = generate_mesh3d(h);
n_nodes = size(nodes, 1);

f_value = @(x) 3*a*pi^2*sin(pi*x(1))*sin(pi*x(2))*sin(pi*x(3)) ...
    - pi*(b(1)*cos(pi*x(1))*sin(pi*x(2))*sin(pi*x(3)) ...
        + b(2)*sin(pi*x(1))*cos(pi*x(2))*sin(pi*x(3)) ...
        + b(3)*sin(pi*x(1))*sin(pi*x(2))*cos(pi*x(3)));

[u_h, ~] = ibeam3d(a, b, f_value, u_true, elements, nodes, n_nodes, nquad, false);

u_ex = u_true(nodes(:,1), nodes(:,2), nodes(:,3));
u_err = abs(u_h - u_ex);

%% Slice plots at x=0.5, y=0.5, z=0.5
slice_dirs = {'x', 'y', 'z'};
slice_vals = {0.5, 0.5, 0.5};
tol = h * 0.6;  % tolerance to grab nodes near slice

for s = 1:3
    figure('Name', sprintf('Solution slice %s=0.5', slice_dirs{s}));

    switch s
        case 1; mask = abs(nodes(:,1) - 0.5) < tol; c1 = 2; c2 = 3; xl = 'y'; yl = 'z';
        case 2; mask = abs(nodes(:,2) - 0.5) < tol; c1 = 1; c2 = 3; xl = 'x'; yl = 'z';
        case 3; mask = abs(nodes(:,3) - 0.5) < tol; c1 = 1; c2 = 2; xl = 'x'; yl = 'y';
    end

    idx  = find(mask);
    xp   = nodes(idx, c1);
    yp   = nodes(idx, c2);
    u_p  = u_h(idx);
    u_ep = u_ex(idx);

    subplot(1,2,1);
    scatter(xp, yp, 40, u_p, 'filled');
    colorbar; axis square; grid on;
    xlabel(xl); ylabel(yl);
    title(sprintf('u_h at %s=0.5', slice_dirs{s}));

    subplot(1,2,2);
    scatter(xp, yp, 40, u_ep, 'filled');
    colorbar; axis square; grid on;
    xlabel(xl); ylabel(yl);
    title(sprintf('u_{true} at %s=0.5', slice_dirs{s}));
end

%% Volume scatter plot colored by error
figure('Name', 'Pointwise error');
scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 20, u_err, 'filled');
colorbar; axis equal; grid on;
xlabel('x'); ylabel('y'); zlabel('z');
title(sprintf('|u_h - u_{true}|, max=%.2e', max(u_err)));
colormap(hot);

%% Isosurface of numerical solution
figure('Name', 'Isosurface u_h = 0.5');
n_plot = k+1;
xg = linspace(0,1,n_plot);
[Xg,Yg,Zg] = meshgrid(xg,xg,xg);

% Interpolate u_h onto regular grid using griddata
U_grid = griddata(nodes(:,1), nodes(:,2), nodes(:,3), u_h, Xg, Yg, Zg);

isoval = 0.5 * max(u_h);
p = patch(isosurface(Xg, Yg, Zg, U_grid, isoval));
isonormals(Xg, Yg, Zg, U_grid, p);
p.FaceColor = [0.2 0.6 1.0];
p.EdgeColor = 'none';
camlight; lighting gouraud;
axis equal; grid on;
xlabel('x'); ylabel('y'); zlabel('z');
title(sprintf('Isosurface u_h = %.2f', isoval));