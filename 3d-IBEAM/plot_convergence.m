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