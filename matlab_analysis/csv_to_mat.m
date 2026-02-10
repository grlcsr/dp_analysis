%% CSV -> MAT bulk converter
% Writes directly in the <bin> directory as:
%   0H.mat
%   <H>_<bit>.mat
clear; clc;

% --------- CONFIG ---------
% We did opendp and ibm with opensll and rap2 with eps 0.1
libraries = {'opendp_mod','opendp', 'ibm'};
bin_names = {'rap2', 'openssl', 'xor', 'mgb_bad'};
mechanism = {'l', 'g'};
epsilons  = 0.1;

% manipulated entropy settings
H_levels  = {"M4", "M8", "M16", "M32", "M64"};
bits_sets = {"clear", "set", "flip", "alternate", "correlation"};

iterations = 1000000;
csv_root  = "D:/";   % input csv dir
out_root  = "./data_mats";   % output dir
% --------------------------

for eps = epsilons
    eps_folder = sprintf("eps_%0.2f", eps);   % e.g., eps_0.01

    for L = 1:numel(libraries)
        library = string(libraries{L});

        % per-library total rows 
        switch lower(library)
            case 'opendp', tot_iterations = 100e6;
            case 'opendp_mod', tot_iterations = 100e6;
            case 'ibm',    tot_iterations = 91e6;
            otherwise,     tot_iterations = [];   % auto-detect if unknown
        end

        for B = 1:numel(bin_names)
            bin = string(bin_names{B});

            % Output BIN directory (flattened target)
            bin_out_dir = fullfile(out_root, eps_folder, library, bin);
            if ~isfolder(bin_out_dir), mkdir(bin_out_dir); end

            % --------- CASE A: baseline 0H ---------
            in_dir_0H = fullfile(csv_root, eps_folder, library + "_csv", bin, "0H");
            for D = 1:numel(mechanism)
                m = string(mechanism{D});
                csv_name = fullfile(in_dir_0H, m + "_eps_" + sprintf('%.3f', eps) + ".csv");
                out_mat  = fullfile(bin_out_dir, "0H.mat");             % <--- NEW NAME/LOCATION

                convert_csv_to_mat_cols23(csv_name, bin_out_dir, out_mat, ...
                    iterations, tot_iterations, eps, library, bin, m, ...
                    "cur_h","", "cur_b","", "overwrite", true);
            end

            % --------- CASE B: manipulated entropy H/bit ---------
            for h = 1:numel(H_levels)
                cur_h = string(H_levels{h});
                for bb = 1:numel(bits_sets)
                    cur_b = string(bits_sets{bb});

                    in_dir_HB = fullfile(csv_root, eps_folder, library + "_csv", bin, cur_h, cur_b);
                    for D = 1:numel(mechanism)
                        m = string(mechanism{D});
                        csv_name = fullfile(in_dir_HB, m + "_eps_" + sprintf('%.3f', eps) + ".csv");
                        out_mat  = fullfile(bin_out_dir, cur_h + "_" + cur_b + ".mat");

                        tot_iterations = iterations;
                        convert_csv_to_mat_cols23(csv_name, bin_out_dir, out_mat, ...
                            iterations, tot_iterations, eps, library, bin, m, ...
                            "cur_h",cur_h, "cur_b",cur_b, "overwrite", true);
                    end
                end
            end
        end
    end
end

fprintf("\nAll conversions finished. Output under: %s\n", out_root);


%% --------- Helper: converter (streams; reads cols #2 & #3; v7.3; creates folders) ---------
function convert_csv_to_mat_cols23(csv_name, out_dir, out_mat, iterations, tot_iterations, eps, library, bin, db, varargin)
    % Optional args
    p = inputParser; p.KeepUnmatched = true;
    addParameter(p, "cur_h", "", @(s)isstring(s) || ischar(s));
    addParameter(p, "cur_b", "", @(s)isstring(s) || ischar(s));
    addParameter(p, "overwrite", true, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    cur_h = string(p.Results.cur_h);
    cur_b = string(p.Results.cur_b);
    overwrite = p.Results.overwrite;

    if ~isfile(csv_name)
        warning("Missing CSV: %s (skipping)", csv_name);
        return;
    end

    if ~isfolder(out_dir), mkdir(out_dir); end
    fprintf("Converting -> %s\n", out_mat);

    % Datastore: header on line 1; select by position (#2 & #3)
    ds = tabularTextDatastore(csv_name, ...
        'ReadVariableNames', true, ...
        'TextType', 'string', ...
        'VariableNamingRule','modify');

    vnames = ds.VariableNames;
    if numel(vnames) < 3
        warning("File has fewer than 3 columns: %s (skipping)", csv_name);
        return;
    end
    ds.SelectedVariableNames = vnames([2 3]);   % always use cols 2 & 3
    ds.ReadSize = iterations;

    % Determine total rows
    if isempty(tot_iterations)
        reset(ds); totN = 0;
        while hasdata(ds); T = read(ds); totN = totN + height(T); end
        reset(ds);
    else
        totN = tot_iterations;
    end

    % Prepare v7.3 MAT and preallocate
    if overwrite && isfile(out_mat), delete(out_mat); end
    m = matfile(out_mat, 'Writable', true);
    m.dp1 = zeros(totN, 1, 'double');   % column 2
    m.dp2 = zeros(totN, 1, 'double');   % column 3

    % metadata
    m.epsilon    = eps;
    m.library    = char(library);
    m.bin        = char(bin);
    m.db_name    = char(db);
    m.iterations = iterations;
    if strlength(cur_h) > 0, m.H_level = char(cur_h); end
    if strlength(cur_b) > 0, m.bit_set = char(cur_b); end

    % Stream & fill
    idx = 1; reset(ds);
    while hasdata(ds)
        T = read(ds);
        n = height(T);
        vv = T.Properties.VariableNames;   % the two selected columns
        v1 = double(T.(vv{1}));
        v2 = double(T.(vv{2}));
        m.dp1(idx:idx+n-1, 1) = v1;
        m.dp2(idx:idx+n-1, 1) = v2;
        idx = idx + n;
    end

    if idx-1 ~= totN
        warning("Row count mismatch for %s (filled %d, expected %d)", out_mat, idx-1, totN);
    end

    fprintf("  Done: %s (rows: %d)\n", out_mat, idx-1);
end
