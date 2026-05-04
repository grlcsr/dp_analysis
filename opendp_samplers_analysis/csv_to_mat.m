% Convert the CSVs from the cks20 sampler runs into .mat files for the
% statistical and plrv analysis

clear; clc;

csv_root  = "./csvs";
out_root  = "./data_mats";
overwrite = true;

% Columns whose values can exceed 2^53 and must stay as strings so they
% don't get silently rounded to double
big_int_columns = ["unif"];

if ~isfolder(out_root), mkdir(out_root); end

trace_cols = ["unif","unif_below","bernoulli","bernoulli_rat", ...
              "bernoulli_exp1","bernoulli_exp","geom_slow","geom_fast","dlap"];

trace_files = {
    "trace_original", trace_cols;
    "trace_taylor",   trace_cols;
    "trace_fixed",    trace_cols;
};

perf_files = {
    "dlap_original", ["d1","d2"];
    "dlap_taylor",   ["d1","d2"];
    "dlap_fixed",    ["d1","d2"];
    "performance",   ["original_ns","taylor_ns","fixed_ns"];
};

traces_runs = list_run_dirs(fullfile(csv_root, "traces"));
perf_runs   = list_run_dirs(fullfile(csv_root, "perf"));

fprintf("Found %d trace run(s) and %d perf run(s)\n\n", ...
        numel(traces_runs), numel(perf_runs));

for r = 1:numel(traces_runs)
    run_name = traces_runs(r);
    src_dir  = fullfile(csv_root, "traces", run_name);
    dst_dir  = fullfile(out_root, "traces", run_name);
    fprintf("=== traces / %s ===\n", run_name);
    convert_directory(src_dir, dst_dir, trace_files, big_int_columns, overwrite);
end

for r = 1:numel(perf_runs)
    run_name = perf_runs(r);
    src_dir  = fullfile(csv_root, "perf", run_name);
    dst_dir  = fullfile(out_root, "perf", run_name);
    fprintf("=== perf / %s ===\n", run_name);
    convert_directory(src_dir, dst_dir, perf_files, big_int_columns, overwrite);
end

fprintf("\nDone. Output under: %s\n", out_root);


function names = list_run_dirs(parent)
    names = strings(0,1);
    if ~isfolder(parent), return; end
    d = dir(parent);
    for i = 1:numel(d)
        if d(i).isdir && startsWith(d(i).name, "run_")
            names(end+1,1) = string(d(i).name);
        end
    end
    names = sort(names);
end

function convert_directory(src_dir, dst_dir, files, big_int_columns, overwrite)
    if ~isfolder(dst_dir), mkdir(dst_dir); end
    for k = 1:size(files, 1)
        base      = string(files{k, 1});
        var_names = files{k, 2};
        csv_path  = fullfile(src_dir, base + ".csv");
        mat_path  = fullfile(dst_dir, base + ".mat");
        convert_csv_to_mat(csv_path, mat_path, var_names, big_int_columns, overwrite);
    end
end

function convert_csv_to_mat(csv_path, mat_path, var_names, big_int_columns, overwrite)
    if ~isfile(csv_path)
        warning("Missing CSV: %s (skipping)", csv_path);
        return;
    end
    out_dir = fileparts(mat_path);
    if strlength(out_dir) > 0 && ~isfolder(out_dir), mkdir(out_dir); end
    fprintf("  Converting %s -> %s\n", csv_path, mat_path);

    if overwrite && isfile(mat_path), delete(mat_path); end

    % Force big-int columns to string before readtable parses them, otherwise
    % they get coerced to double and lose precision
    opts = detectImportOptions(csv_path, ...
        'VariableNamingRule', 'preserve', ...
        'TextType', 'string');
    for i = 1:numel(big_int_columns)
        name = big_int_columns(i);
        if any(strcmp(opts.VariableNames, name))
            opts = setvartype(opts, name, 'string');
        end
    end

    T = readtable(csv_path, opts);

    n_cols = numel(var_names);
    if width(T) < n_cols
        warning("File has fewer than %d columns: %s (skipping)", n_cols, csv_path);
        return;
    end

    s = struct();
    for c = 1:n_cols
        col = T.(T.Properties.VariableNames{c});
        s.(char(var_names(c))) = parse_column(col);
    end

    s.source_csv = char(csv_path);
    s.n_rows     = height(T);
    s.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

    save(mat_path, '-struct', 's', '-v7.3');
    fprintf("    Done: %s (rows: %d, cols: %d)\n", mat_path, height(T), n_cols);
end

function out = parse_column(col)
    if islogical(col)
        out = double(col);
        return;
    end
    if isnumeric(col)
        out = double(col);
        return;
    end

    s = string(col);

    % true/false columns map to 0/1
    is_true = (s == "true");
    if all(is_true | s == "false" | s == "")
        out = double(is_true);
        return;
    end

    % > 17 digits can't have round-tripped through a double without losing
    % precision, so keep these as strings
    if any(strlength(s) > 17)
        out = s;
        return;
    end

    out = str2double(s);
end