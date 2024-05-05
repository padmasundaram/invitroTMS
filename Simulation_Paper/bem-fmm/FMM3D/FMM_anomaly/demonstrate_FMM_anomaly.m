%% Demonstrate anomaly in FMM3d for Laplacian interaction library
disp("Demonstrate anomaly in FMM3d for Laplacian interaction library.");
disp("Run this script once on Windows using the precompiled lfmm3d mexa shipped with BEM-FMM (included here) and once on a Linux-based OS using a fresh compilation (without FAST_KER).");
disp("Electric field and potential may be computed directly on machines with enough working memory. Use the flag 'is_direct' in the script. Otherwise, a precomputed exact solution will be loaded." + newline);

is_direct = false;

%% Load input
load("FMM_anomaly_input", 'srcinfo', 'pg');
prec = 1e-2;

%% Call FMM3d library
U = lfmm3d(prec, srcinfo, pg);

%% Compute electric field and potential directly
if is_direct
    disp("Compute electric potential and field directly." + newline);
    
    P = srcinfo.sources';
    C = srcinfo.charges';

    % Electric potential
    diff_x = P(:, 1) - P(:, 1)';                    % r_x-r'_x
    diff_y = P(:, 2) - P(:, 2)';                    % r_y-r'_y
    diff_z = P(:, 3) - P(:, 3)';                    % r_z-r'_z
    dist = sqrt(diff_x.^2 + diff_y.^2 + diff_z.^2); % |r-r'|
    ker = 1./dist;                                  % 1/|r-r'|
    ker(1:size(P, 1)+1:end) = 0;                    % Set r'=r terms to 0
    potential_direct = (ker*C);                     % Compute sum_i c_i/|r-r_i'|
    
    % Electric field
    ker_x = diff_x.*ker.^3;                         % r_x-r'_x/|r-r'|^3
    ker_x(1:size(P, 1)+1:end) = 0;                  % Set r'=r terms to 0
    ker_y = diff_y.*ker.^3;                         % r_y-r'_y/|r-r'|^3
    ker_y(1:size(P, 1)+1:end) = 0;                  % Set r'=r terms to 0
    ker_z = diff_z.*ker.^3;                         % r_z-r'_z/|r-r'|^3
    ker_z(1:size(P, 1)+1:end) = 0;                  % Set r'=r terms to 0
    field_direct = zeros(size(P, 1), 3);
    field_direct(:, 1) = ker_x*C;                   % Compute sum_i c_i (r_x-r'_x)/|r-r_i'|^3
    field_direct(:, 2) = ker_y*C;                   % Compute sum_i c_i (r_y-r'_y)/|r-r_i'|^3
    field_direct(:, 3) = ker_z*C;                   % Compute sum_i c_i (r_z-r'_z)/|r-r_i'|^3

    save("FMM_anomaly_exact_precomputed.mat", 'potential_direct', 'field_direct');
else    
    load("FMM_anomaly_exact_precomputed.mat");    
end

disp("Compare FMM and direct results.");
disp("Euclidean distance between the two potential vectors:");
disp(norm(potential_direct - U.pot'));
disp("Euclidean distance between the two field matrices (computed over all entries):");
disp(sqrt(sum((field_direct - (-U.grad')).^2, 'all')));
disp("Euclidean distance between the two potential vectors normalized by norm of directly computed potential:");
disp(norm(potential_direct - U.pot')./norm(potential_direct));
disp("Euclidean distance between the two field matrices (computed over all entries) normalized by norm (over all entries) of directly computed field:");
disp(sqrt(sum((field_direct - (-U.grad')).^2, 'all'))./sqrt(sum(field_direct.^2, 'all')));
