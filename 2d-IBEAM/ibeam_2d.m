addpath('PolyMesher')
close all;

b = [0; 1];
k_values = 3:7;
a_values = (1/10).^(1:4);

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
        %f_value = @(x,y) 0 + a*2*sin(x+y) - (b(1) + b(2)) * cos(x + y);
        f_value = @(x,y) 2*a*pi^2*sin(pi*x).*sin(pi*y) - ...
            pi*(b(1)*cos(pi*x).*sin(pi*y) + b(2)*sin(pi*x).*cos(pi*y));
        u_true = @(x,y) 0 + sin(pi*x)*sin(pi*y);
        u_x = @(x,y) 0 + pi*cos(pi*x)*sin(pi*y);
        u_y = @(x,y) 0 + pi*cos(pi*y)*sin(pi*x);

        [u_ibeam, A] = ibeam(a, b, f_value, u_true, elements, nodes, n_nodes, nquad);

        % [X,Y] = meshgrid(linspace(0,1,n));
        % figure('Name', sprintf('Graph of Solution: n=%d, alpha=%.3f', n, a))
        % surf(X,Y, reshape(u_ibeam, n, n)); hold on
        % xlabel('x'); ylabel('y'); zlabel('u(x,y)');
        % axis square
        % grid on
        % title(sprintf('Solution for n = %d, \\alpha = %.3f', n, a))

        H1_ibeam_k(idx,a_idx) = H1_semi_error(u_ibeam, u_x, u_y, elements, nodes, nquad);
        L2_ibeam_k(idx,a_idx) = L2_error(u_ibeam, u_true, elements, nodes, nquad);
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

figure('Name', 'L^2 norm Error')
loglog(h_values, L2_ibeam_k(:,1), 'LineWidth', 1.5); hold on
loglog(h_values, L2_ibeam_k(:,2), 'LineWidth', 1.5)
loglog(h_values, L2_ibeam_k(:,3), 'LineWidth', 1.5)
loglog(h_values, L2_ibeam_k(:,4), 'LineWidth', 1.5)
loglog(h_values, C1*h_values, '--k', 'LineWidth', 1.2)
loglog(h_values, C2*h_values.^2, '--r', 'LineWidth', 1.2)
xlabel('h'); ylabel('error'); title('L^2 Norm Error')
legend('\alpha = 0.1', '\alpha = 0.01', '\alpha = 1E-3', '\alpha = 1E-4', 'O(h)', 'O(h^2)', 'Location', 'best'); grid on

C1 = H1_ibeam_k(1,1) / h_values(1);
C2 = H1_ibeam_k(1,1) / h_values(1)^2;

figure('Name', 'H^1 Semi-Norm Error')
loglog(h_values, H1_ibeam_k(:,1), 'LineWidth', 1.5); hold on
loglog(h_values, H1_ibeam_k(:,2), 'LineWidth', 1.5)
loglog(h_values, H1_ibeam_k(:,3), 'LineWidth', 1.5)
loglog(h_values, H1_ibeam_k(:,4), 'LineWidth', 1.5)
loglog(h_values, C1*h_values, '--k', 'LineWidth', 1.2)
loglog(h_values, C2*h_values.^2, '--r', 'LineWidth', 1.2)
xlabel('h'); ylabel('error'); title('H^1 Semi-norm Error')
legend('\alpha = 0.1', '\alpha = 0.01', '\alpha = 1E-3', '\alpha = 1E-4', 'O(h)', 'O(h^2)', 'Location', 'best'); grid on
