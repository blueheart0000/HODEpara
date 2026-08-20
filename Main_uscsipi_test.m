%% Benchmark: Multi-Threshold Segmentation on test_images/
% Evaluates all algorithms on all .mat images in the test_images/
% directory. Each .mat file must contain a variable 'img' (uint8
% grayscale matrix). Image name is derived from the filename.
%
% Usage:  >> benchmark_uscsipi
% Output: console summary + figs1/benchmark_uscsipi_results.txt

clear; clc; close all;
addpath(pwd);

rng(1400200, 'twister');

%% ======================== Configuration ========================
K         = 4;           % number of thresholds
L         = 256;         % gray levels
method    = 'otsu';      % 'otsu' or 'kapur'
num_runs  = 25;          % runs per image per algorithm
pop_size  = 50;
max_fes   = 140000;

%% ======================== Test Set ========================
% Dynamically loads all .mat files from test_images/ directory.
% Each .mat file must contain a variable 'img' (uint8 grayscale matrix).
% The image name is derived from the filename (without extension).

test_dir = fullfile(pwd, 'test_images');
if ~exist(test_dir, 'dir'), mkdir(test_dir); end

% Scan all .mat files in test_images/
mat_files = dir(fullfile(test_dir, '*.mat'));
num_images = length(mat_files);
if num_images == 0
    error('No .mat files found in %s.\n  Please place your test images as .mat files with variable ''img''.', test_dir);
end

fprintf('Loading %d test images from %s ...\n', num_images, test_dir);
test_set = cell(num_images, 2);

for im = 1:num_images
    fname = mat_files(im).name;
    [~, name, ~] = fileparts(fname);
    img_path = fullfile(test_dir, fname);

    s = load(img_path, 'img');
    img_uint8 = s.img;
    test_set{im, 1} = name;
    test_set{im, 2} = @(N) img_uint8;
    fprintf('  %d. %s  (%d x %d)\n', im, name, size(img_uint8,1), size(img_uint8,2));
end
fprintf('Loaded %d test images.\n\n', num_images);

%% ======================== Algorithms ========================
algorithms = {
    'DE\_rand\_1',  @(obj_func, dim, X0) run_de_rand1(obj_func, dim, pop_size, max_fes, X0);
    'SHADE\_call',  @(obj_func, dim, X0) run_shade_call(obj_func, dim, pop_size, max_fes, X0);
    'L\_SHADE',     @(obj_func, dim, X0) run_lshade(obj_func, dim, pop_size, max_fes);
    'jSO',          @(obj_func, dim, X0) run_jso(obj_func, dim, pop_size, max_fes, X0);
    'EhdDE',        @(obj_func, dim, X0) run_ehdde(obj_func, dim, pop_size, max_fes, X0);
    'HODE',         @(obj_func, dim, X0) run_hode(obj_func, dim, pop_size, max_fes, X0);
};
num_algs = size(algorithms, 1);

fprintf('================================================================\n');
fprintf('  USC-SIPI Multi-Threshold Segmentation Benchmark\n');
fprintf('  Test set: %d images, K=%d, Method=%s\n', num_images, K, method);
fprintf('  Algorithms: %d, Runs/img: %d, MaxFEs: %d\n', num_algs, num_runs, max_fes);
fprintf('================================================================\n\n');

%% ======================== Storage ========================
% [num_images x num_algs] with mean and std
img_fit_mean = zeros(num_images, num_algs);
img_fit_std  = zeros(num_images, num_algs);
img_psnr_mean = zeros(num_images, num_algs);
img_psnr_std  = zeros(num_images, num_algs);
img_ssim_mean = zeros(num_images, num_algs);
img_ssim_std  = zeros(num_images, num_algs);
img_names     = cell(num_images, 1);
img_sizes     = zeros(num_images, 2);

%% ======================== Run on Each Image ========================
out_dir = 'figs1';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

for im = 1:num_images
    img_name = test_set{im, 1};
    img_fn   = test_set{im, 2};
    img_uint8 = img_fn(256);
    img_names{im} = img_name;
    img_sizes(im,:) = size(img_uint8);

    hist_counts = imhist(img_uint8, L)';
    dim = K;

    switch lower(method)
        case 'otsu',  obj_func = @(x) otsu_multi_thresh(x, hist_counts, L);  obj_name = 'Otsu';
        case 'kapur', obj_func = @(x) kapur_multi_thresh(x, hist_counts, L); obj_name = 'Kapur';
    end

    % Per-image per-algorithm raw data: [num_runs x num_algs]
    img_fit_raw = zeros(num_runs, num_algs);
    img_ps_raw  = zeros(num_runs, num_algs);
    img_ss_raw  = zeros(num_runs, num_algs);
    img_thresh_raw = cell(num_runs, num_algs);

    fprintf('--- Image %d/%d: %s (%d x %d) ---\n', im, num_images, img_name, size(img_uint8,1), size(img_uint8,2));

    parfor run = 1:num_runs
        X0 = rand(pop_size, dim);  % shared initial population for all algorithms on this run
        tmp_fit = zeros(1, num_algs);
        tmp_ps  = zeros(1, num_algs);
        tmp_ss  = zeros(1, num_algs);
        tmp_thresh = cell(1, num_algs);
        for a = 1:num_algs
            alg_func = algorithms{a, 2};
            [best_sol, best_fit] = alg_func(obj_func, dim, X0);

            thresholds = sort(max(1, min(L-2, round(best_sol(:)' * (L-2)) + 1)));
            % Guard: if thresholds empty or all same, use evenly spaced defaults
            if isempty(thresholds) || length(unique(thresholds)) < K
                thresholds = round(linspace(1, L-2, K));
            end
            seg = apply_thresholds(img_uint8, thresholds, L);

            tmp_fit(a) = best_fit;
            tmp_thresh{a} = thresholds;

            mse_val = mean((double(seg(:)) - double(img_uint8(:))).^2);
            if mse_val > 0
                tmp_ps(a) = 10 * log10(255^2 / mse_val);
            end
            tmp_ss(a) = ssim_index_local(double(seg), double(img_uint8));
        end
        img_fit_raw(run, :) = tmp_fit;
        img_ps_raw(run, :)  = tmp_ps;
        img_ss_raw(run, :)  = tmp_ss;
        img_thresh_raw(run, :) = tmp_thresh;
    end

    % Print per-algorithm statistics for this image
    for a = 1:num_algs
        fprintf('  %s: fit=%.4e+-%.4e, PSNR=%.2f+-%.2f\n', ...
            algorithms{a,1}, ...
            abs(mean(img_fit_raw(:,a))), std(img_fit_raw(:,a)), ...
            mean(img_ps_raw(:,a)), std(img_ps_raw(:,a)));
    end

    % Store mean & std
    img_fit_mean(im,:)  = mean(img_fit_raw, 1);
    img_fit_std(im,:)   = std(img_fit_raw, 0, 1);
    img_psnr_mean(im,:) = mean(img_ps_raw, 1);
    img_psnr_std(im,:)  = std(img_ps_raw, 0, 1);
    img_ssim_mean(im,:) = mean(img_ss_raw, 1);
    img_ssim_std(im,:)  = std(img_ss_raw, 0, 1);

    % Save segmentation images for this image
    img_out_dir = fullfile(out_dir, img_name);
    if ~exist(img_out_dir, 'dir'), mkdir(img_out_dir); end

    % Save original image
    imwrite(img_uint8, fullfile(img_out_dir, 'original.png'));

    % For each algorithm, pick the run with median fitness and save segmentation
    for a = 1:num_algs
        fit_vals = img_fit_raw(:, a);
        [~, sort_idx] = sort(fit_vals);
        median_run = sort_idx(round(num_runs/2));
        best_thresh = img_thresh_raw{median_run, a};
        % Guard against empty/collapsed thresholds
        if isempty(best_thresh) || length(unique(best_thresh)) < K
            fprintf('  WARNING: %s bad thresholds = [%s], using defaults\n', ...
                strrep(algorithms{a,1}, '\_', '_'), sprintf('%d ', best_thresh));
            best_thresh = round(linspace(1, L-2, K));
        end
        seg_img = apply_thresholds(img_uint8, best_thresh, L);
        alg_name_short = strrep(algorithms{a,1}, '\_', '_');
        imwrite(seg_img, fullfile(img_out_dir, [alg_name_short, '.png']));
    end
    fprintf('  Segmentation images saved to %s\n', img_out_dir);
end

%% ======================== Print Results Table ========================
fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  RESULTS: Fitness (mean +- std) for each image\n');
fprintf('  (smaller fitness = better objective value)\n');
fprintf('================================================================\n\n');

abs_fit_mean = abs(img_fit_mean);
abs_fit_std  = img_fit_std;

for im = 1:num_images
    fprintf('--- %s (%d x %d) ---\n', img_names{im}, img_sizes(im,1), img_sizes(im,2));
    fprintf('%-12s', 'Algorithm');
    fprintf(' %-24s', 'Fitness');
    fprintf(' %-16s', 'PSNR(dB)');
    fprintf(' %-12s\n', 'SSIM');
    fprintf('%s\n', repmat('-', 12+24+16+12, 1));
    for a = 1:num_algs
        fprintf('%-12s', algorithms{a,1});
        fprintf(' %-10.4e +- %-8.4e', abs_fit_mean(im,a), abs_fit_std(im,a));
        fprintf(' %-8.2f +- %-5.2f', img_psnr_mean(im,a), img_psnr_std(im,a));
        fprintf(' %-6.4f +- %-5.4f\n', img_ssim_mean(im,a), img_ssim_std(im,a));
    end
    fprintf('\n');
end

%% ======================== Friedman Ranking ========================
fprintf('================================================================\n');
fprintf('  FRIEDMAN TEST (overall ranking across all %d images)\n', num_images);
fprintf('================================================================\n\n');

[p_f, ~, stats_f] = friedman(abs_fit_mean, 1, 'off');
ranks = stats_f.meanranks;
[~, idx] = sort(ranks);

fprintf('  Friedman p-value = %.6f\n', p_f);
fprintf('  Ranking based on Fitness (lower rank = better):\n');
fprintf('  %-12s %-12s\n', 'Algorithm', 'Mean Rank');
fprintf('  %s\n', repmat('-', 25, 1));
for r = 1:num_algs
    fprintf('  %d. %-10s %.3f\n', r, algorithms{idx(r),1}, ranks(idx(r)));
end

%% ======================== Save Results to Excel ========================
xlsx_file = fullfile(out_dir, 'benchmark_uscsipi_results.xlsx');

% Delete old file if exists
if exist(xlsx_file, 'file'), delete(xlsx_file); end

% Column arrangement:
%   Image  |  Fitness_mean (DE .. HODE)  |  Fitness_std (DE .. HODE)
%   |  PSNR_mean (DE .. HODE)  |  PSNR_std (DE .. HODE)
%   |  SSIM_mean (DE .. HODE)  |  SSIM_std (DE .. HODE)
metrics = {'Fitness_mean','Fitness_std','PSNR_mean','PSNR_std','SSIM_mean','SSIM_std'};
alg_clean = cell(num_algs, 1);
for a = 1:num_algs
    alg_clean{a} = strrep(algorithms{a,1}, '\_', '_');
end

header = {'Image'};
for m = 1:length(metrics)
    for a = 1:num_algs
        header{end+1} = [alg_clean{a} '_' metrics{m}];
    end
end

% Build data: one row per image
data_cell = cell(num_images+1, length(header));
for c = 1:length(header)
    data_cell{1,c} = header{c};
end
for im = 1:num_images
    data_cell{im+1, 1} = img_names{im};
    col = 2;
    % Fitness_mean
    for a = 1:num_algs
        data_cell{im+1, col} = abs_fit_mean(im,a); col = col + 1;
    end
    % Fitness_std
    for a = 1:num_algs
        data_cell{im+1, col} = abs_fit_std(im,a); col = col + 1;
    end
    % PSNR_mean
    for a = 1:num_algs
        data_cell{im+1, col} = img_psnr_mean(im,a); col = col + 1;
    end
    % PSNR_std
    for a = 1:num_algs
        data_cell{im+1, col} = img_psnr_std(im,a); col = col + 1;
    end
    % SSIM_mean
    for a = 1:num_algs
        data_cell{im+1, col} = img_ssim_mean(im,a); col = col + 1;
    end
    % SSIM_std
    for a = 1:num_algs
        data_cell{im+1, col} = img_ssim_std(im,a); col = col + 1;
    end
end
xlswrite(xlsx_file, data_cell, 'All_Results');

% Summary sheet: Friedman rankings
rank_cell = cell(num_algs+1, 2);
rank_cell{1,1} = 'Algorithm'; rank_cell{1,2} = 'Friedman_Rank';
for a = 1:num_algs
    rank_cell{a+1,1} = algorithms{a,1};
    rank_cell{a+1,2} = ranks(a);
end
xlswrite(xlsx_file, rank_cell, 'Friedman_Rank');

fprintf('\nExcel results saved to: %s\n', xlsx_file);

%% ======================== Save Results to Text (backup) ========================
txt_file = fullfile(out_dir, 'benchmark_uscsipi_results.txt');
fid = fopen(txt_file, 'w');
fprintf(fid, 'USC-SIPI Benchmark Results\n');
fprintf(fid, 'K=%d, Method=%s, Runs=%d, MaxFEs=%d\n\n', K, method, num_runs, max_fes);
for im = 1:num_images
    fprintf(fid, '--- %s (%d x %d) ---\n', img_names{im}, img_sizes(im,1), img_sizes(im,2));
    fprintf(fid, '%-12s %-24s %-16s %-12s\n', 'Algorithm', 'Fitness', 'PSNR(dB)', 'SSIM');
    for a = 1:num_algs
        fprintf(fid, '%-12s %-10.4e+-%-10.4e %-8.2f+-%-6.2f %-6.4f+-%-6.4f\n', ...
            algorithms{a,1}, abs_fit_mean(im,a), abs_fit_std(im,a), ...
            img_psnr_mean(im,a), img_psnr_std(im,a), ...
            img_ssim_mean(im,a), img_ssim_std(im,a));
    end
    fprintf(fid, '\n');
end
fprintf(fid, 'Friedman p-value = %.6f\n', p_f);
fprintf(fid, 'Algorithm Rankings:\n');
for r = 1:num_algs
    fprintf(fid, '  %d. %s (rank %.3f)\n', r, algorithms{idx(r),1}, ranks(idx(r)));
end
fclose(fid);
fprintf('Text results saved to: %s\n', txt_file);

%% ======================== Visualization ========================
figure('Name', 'USC-SIPI Benchmark', 'Position', [50 50 1400 900]);

% Subplot 1: Friedman rankings
subplot(2,3,1);
bar(ranks);
set(gca, 'XTickLabel', algorithms(:,1));
ylabel('Mean Rank');
title('Friedman Rank (lower=better)');
grid on;

% Subplot 2: Average PSNR across all images
subplot(2,3,2);
mean_psnr = mean(img_psnr_mean, 1);
std_psnr  = mean(img_psnr_std, 1);
barwitherr(std_psnr, mean_psnr);
set(gca, 'XTickLabel', algorithms(:,1));
ylabel('PSNR (dB)');
title('Mean PSNR Across All Images');
grid on;

% Subplot 3: Average SSIM
subplot(2,3,3);
mean_ssim = mean(img_ssim_mean, 1);
std_ssim  = mean(img_ssim_std, 1);
barwitherr(std_ssim, mean_ssim);
set(gca, 'XTickLabel', algorithms(:,1));
ylabel('SSIM');
title('Mean SSIM Across All Images');
grid on;

% Subplot 4: Fitness heatmap
subplot(2,3,4);
imagesc(abs_fit_mean);
colorbar;
xlabel('Algorithm'); ylabel('Image');
set(gca, 'XTick', 1:num_algs, 'XTickLabel', algorithms(:,1));
set(gca, 'YTick', 1:num_images, 'YTickLabel', img_names);
title('Fitness Heatmap');

% Subplot 5: PSNR heatmap
subplot(2,3,5);
imagesc(img_psnr_mean);
colorbar;
xlabel('Algorithm'); ylabel('Image');
set(gca, 'XTick', 1:num_algs, 'XTickLabel', algorithms(:,1));
set(gca, 'YTick', 1:num_images, 'YTickLabel', img_names);
title('PSNR Heatmap (dB)');

% Subplot 6: Average fitness bar
subplot(2,3,6);
mean_fit_all = mean(abs_fit_mean, 1);
std_fit_all  = mean(abs_fit_std, 1);
barwitherr(std_fit_all, mean_fit_all);
set(gca, 'XTickLabel', algorithms(:,1));
ylabel('Mean Fitness');
title('Average Fitness Across All Images');
grid on;

annotation('textbox', [0, 0.95, 1, 0.05], ...
    'String', sprintf('USC-SIPI Benchmark: %d Images, K=%d, %s, %dx%d runs', ...
    num_images, K, method, num_runs, num_algs), ...
    'HorizontalAlignment', 'center', 'FontSize', 13, ...
    'FontWeight', 'bold', 'EdgeColor', 'none');

%% ======================== Helper Functions ========================

function s = apply_thresholds(img, thresholds, L)
    boundaries = [0, thresholds, L-1];
    s = zeros(size(img), 'uint8');
    for k = 1:(length(thresholds)+1)
        lo = boundaries(k); hi = boundaries(k+1);
        mid = round((lo + hi) / 2);
        s((img >= lo) & (img <= hi)) = mid;
    end
end

function mssim = ssim_index_local(img1, img2)
    K1=0.01; K2=0.03; L=255;
    C1=(K1*L)^2; C2=(K2*L)^2;
    mu1=mean(img1(:)); mu2=mean(img2(:));
    s1=var(img1(:)); s2=var(img2(:));
    s12=mean((img1(:)-mu1).*(img2(:)-mu2));
    mssim=((2*mu1*mu2+C1)*(2*s12+C2))/((mu1^2+mu2^2+C1)*(s1+s2+C2));
end

function h = barwitherr(err, vals)
    h = bar(vals);
    hold on;
    x = 1:length(vals);
    errorbar(x, vals, err, 'k', 'LineStyle', 'none', 'LineWidth', 1);
    hold off;
end

%% ======================== Algorithm Wrappers ========================

function [best_sol, best_fit] = run_de_rand1(obj_func, dim, pop_size, max_fes, X0)
    max_iter=ceil(max_fes/(2*pop_size));
    [best_sol, best_fit, ~] = DE_rand_1(X0, max_iter, pop_size, 0, 1, dim, obj_func);
end

function [best_sol, best_fit] = run_shade_call(obj_func, dim, pop_size, max_fes, X0)
    [best_sol, best_fit]=SHADE_call(obj_func, 0, 1, dim, max_fes, pop_size, X0);
end

function [best_sol, best_fit] = run_lshade(obj_func, dim, pop_size, max_fes, ~)
    xl=[zeros(dim,1), ones(dim,1)];
    options=[pop_size, max_fes, 0];
    [best_sol, best_fit]=L_SHADE(xl, options, obj_func);
end

function [best_sol, best_fit] = run_jso(obj_func, dim, pop_size, max_fes, X0)
    [best_sol, best_fit]=jSO(obj_func, 0, 1, dim, max_fes, pop_size, X0);
end

function [best_sol, best_fit] = run_ehdde(obj_func, dim, pop_size, max_fes, X0)
    max_iter=ceil(max_fes/(6*pop_size));
    [best_sol, best_fit, ~]=EhdDE(X0, max_iter, pop_size, 0, 1, dim, obj_func);
end

function [best_sol, best_fit] = run_hode(obj_func, dim, pop_size, max_fes, X0)
    max_iter=ceil(max_fes/(7*pop_size));
    [best_fit, best_sol, ~]=HODE_para_NoEES(obj_func, pop_size, max_iter, 0, 1, dim, X0);
end
