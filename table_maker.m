% =========================================================
%  make_convergence_tables.m
% =========================================================

clear; clc;

% ----------------------------------------------------------
% 1. CONFIGURE
% ----------------------------------------------------------
experiments_2D = {
%    'squares',             'squares.mat';
%    'quadrilaterals',      'quads.mat';
    'gen.~cvx.~polygons',  'polygons.mat';
};

experiments_3D = {
%    'tetrahedrons',  'tets.mat';
%    'cubes',         'cubes.mat';
};

% ----------------------------------------------------------
% 2. Load all data
% ----------------------------------------------------------
function D = load_experiment(label, fname)
    D = struct('label', label, 'fname', fname, ...
               'ok', false, 'h', [], 'L2', [], 'H1', [], 'a', []);
    if ~isfile(fname)
        fprintf('  *** File not found: %s ***\n', fname);
        return;
    end
    S = load(fname);
    required = {'h_values','L2_err','H1_err','a_values'};
    missing  = required(~isfield(S, required));
    if ~isempty(missing)
        fprintf('  *** %s missing variables: %s ***\n', fname, strjoin(missing,', '));
        fprintf('      Found: %s\n', strjoin(fieldnames(S)', ', '));
        return;
    end
    [h_sorted, idx] = sort(S.h_values(:), 'ascend');  % ascending h
    D.h  = h_sorted;
    D.L2 = S.L2_err(idx, :);
    D.H1 = S.H1_err(idx, :);
    D.a  = S.a_values(:)';
    D.ok = true;
end

data_2D = cell(size(experiments_2D, 1), 1);
for i = 1:size(experiments_2D, 1)
    data_2D{i} = load_experiment(experiments_2D{i,1}, experiments_2D{i,2});
end

data_3D = cell(size(experiments_3D, 1), 1);
for i = 1:size(experiments_3D, 1)
    data_3D{i} = load_experiment(experiments_3D{i,1}, experiments_3D{i,2});
end

% Infer a_values from first successful load
a_values = [];
all_data = [data_2D; data_3D];
for i = 1:numel(all_data)
    if all_data{i}.ok
        a_values = all_data{i}.a;
        break;
    end
end
if isempty(a_values)
    error('No data loaded successfully.');
end

% ----------------------------------------------------------
% 3. Helpers
% ----------------------------------------------------------
function orders = conv_orders(h_vec, err_vec)
    orders = nan(size(h_vec));
    for k = 2:numel(h_vec)
        if err_vec(k) > 0 && err_vec(k-1) > 0
            orders(k) = log(err_vec(k-1)/err_vec(k)) / ...
                        log(h_vec(k-1)/h_vec(k));
        end
    end
end

function s = fmt_ord(x)
    if isnan(x), s = '{-}'; else, s = sprintf('%.2f', x); end
end

function s = fmt_err(x)
    x = x(1);   % ensure scalar
    if x == 0, s = '0'; return; end
    exp_n = floor(log10(abs(x)));
    mant  = x ./ 10.^exp_n;
    if abs(mant) >= 9.9995
        mant  = mant ./ 10;
        exp_n = exp_n + 1;
    end
    if exp_n >= 0
        s = sprintf('%.3f\\text{e}+%d', mant, exp_n);
    else
        s = sprintf('%.3f\\text{e}%d', mant, exp_n);
    end
end

% ----------------------------------------------------------
% 4. Print one LaTeX table
% ----------------------------------------------------------
function print_latex_table(fid, dim_label, data_list, a_val, a_idx)
    n_exp = numel(data_list);

    % Column spec
    col_spec = 'l';
    for i = 1:n_exp
        col_spec = [col_spec, ' | cccc'];
    end

    fprintf(fid, '\\begin{table}[h]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{%s mesh refinement, $a = %g$, ', dim_label, a_val);
    fprintf(fid, '$e_h = u_I - u_h$}\n');
    % resizebox scales to textwidth while preserving aspect ratio
    fprintf(fid, '\\resizebox{\\textwidth}{!}{%%\n');
    fprintf(fid, '\\begin{tabular}{%s}\n', col_spec);
    fprintf(fid, '\\hline\n');

    % Header row 1: mesh type labels
    fprintf(fid, '$h$');
    for i = 1:n_exp
        fprintf(fid, ' & \\multicolumn{4}{c}{%s}', data_list{i}.label);
    end
    fprintf(fid, ' \\\\\n');

    % Cmidrules
    for i = 1:n_exp
        sc = 2 + (i-1)*4;
        ec = sc + 3;
        fprintf(fid, '\\cmidrule(lr){%d-%d}', sc, ec);
    end
    fprintf(fid, '\n');

    % Header row 2: sub-column names
    for i = 1:n_exp
        fprintf(fid, ' & $\\|e_h\\|_{L^2}$ & Ord & $|e_h|_{H^1}$ & Ord');
    end
    fprintf(fid, ' \\\\\n\\hline\n');

    % Unified k values (avoids float comparison duplicates)
    all_k = [];
    for i = 1:n_exp
        if data_list{i}.ok
            k_vals = -round(log2(data_list{i}.h));
            all_k  = union(all_k, k_vals);
        end
    end
    all_k = sort(all_k, 'ascend');  % ascending k = coarse to fine

    % Data rows
    for j = 1:numel(all_k)
        kj = all_k(j);
        fprintf(fid, '$2^{-%d}$', kj);
        for i = 1:n_exp
            D = data_list{i};
            if ~D.ok
                fprintf(fid, ' & -- & -- & -- & --');
                continue;
            end
            k_vals = -round(log2(D.h));
            row    = find(k_vals == kj);
            if isempty(row)
                fprintf(fid, ' & ${-}$ & ${-}$ & ${-}$ & ${-}$');
                continue;
            end
            row = row(1);
            l2_ord = conv_orders(D.h, D.L2(:, a_idx));
            h1_ord = conv_orders(D.h, D.H1(:, a_idx));
            fprintf(fid, ' & $%s$ & %s & $%s$ & %s', ...
                fmt_err(D.L2(row, a_idx)), fmt_ord(l2_ord(row)), ...
                fmt_err(D.H1(row, a_idx)), fmt_ord(h1_ord(row)));
        end
        fprintf(fid, ' \\\\\n');
    end

    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '}%%  end resizebox\n');
    fprintf(fid, '\\end{table}\n');
end

% ----------------------------------------------------------
% 5. Write .tex file
% ----------------------------------------------------------
fid = fopen('convergence_tables.tex', 'w');

fprintf(fid, '\\documentclass[a4paper]{article}\n');
fprintf(fid, '\\usepackage[a4paper, margin=2cm]{geometry}\n');
fprintf(fid, '\\usepackage{booktabs}\n');
fprintf(fid, '\\usepackage{amsmath}\n');
fprintf(fid, '\\begin{document}\n\n');

for ai = 1:numel(a_values)
    print_latex_table(fid, '2D', data_2D, a_values(ai), ai);
    fprintf(fid, '\n\\bigskip\n\n');
    print_latex_table(fid, '3D', data_3D, a_values(ai), ai);
    fprintf(fid, '\n\\clearpage\n\n');
end

fprintf(fid, '\\end{document}\n');
fclose(fid);

fprintf('Done! Output written to convergence_tables.tex\n');