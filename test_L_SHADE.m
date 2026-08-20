% test_L_SHADE.m
% Test L-SHADE on all 30 CEC2017 benchmark functions (F1-F30)
clear; clc;
rng('default')


func_list = 1:12;  % CEC2017 F1 to F30

% Result storage
f_mean   = zeros(1, length(func_list));
f_std    = zeros(1, length(func_list));
f_best   = zeros(1, length(func_list));
f_worst  = zeros(1, length(func_list));


f_mean_1  = zeros(1, length(func_list));
f_std_1   = zeros(1, length(func_list));
f_best_1  = zeros(1, length(func_list));
f_worst_1 = zeros(1, length(func_list));


f_mean_2  = zeros(1, length(func_list));
f_std_2   = zeros(1, length(func_list));
f_best_2  = zeros(1, length(func_list));
f_worst_2 = zeros(1, length(func_list));

f_mean_3  = zeros(1, length(func_list));
f_std_3   = zeros(1, length(func_list));
f_best_3  = zeros(1, length(func_list));
f_worst_3 = zeros(1, length(func_list));

f_mean_4  = zeros(1, length(func_list));
f_std_4   = zeros(1, length(func_list));
f_best_4  = zeros(1, length(func_list));
f_worst_4 = zeros(1, length(func_list));

f_mean_5  = zeros(1, length(func_list));
f_std_5   = zeros(1, length(func_list));
f_best_5  = zeros(1, length(func_list));
f_worst_5 = zeros(1, length(func_list));

tic;
N        = 100;%min(100, 18 * dim);
N_runs   = 25;  % number of independent runs
D        = 20;
% FE budget
max_fes  = 1400200;

% Per-generation FE costs for generation-based algorithms:
%   DE_rand_1: 2N/gen
%   EhdDE:     6N/gen
%   HODE:      7N/gen + 2N initial

% Compute iteration counts corresponding to max_fes for each algorithm
T       = round(max_fes / (2 * N));   % DE_rand_1 iterations    → 7001
T_ehdde = round(max_fes / (6 * N));   % EhdDE iterations        → 2334
T_hode  = round((max_fes - 2*N) / (7 * N));  % HODE iterations  → 2000



parfor fi = func_list
    % Get CEC2017 benchmark function info
    [lb, ub, dim, fobj] = Get_CEC2022_Functions_details(strcat('F', num2str(fi)),D);

    % Create initial range: dim x 2 matrix, each row = [lb, ub]
    if isscalar(lb)
        xl = repmat([lb, ub], dim, 1);
    else
        xl = [lb(:), ub(:)];
    end
    
    x0 = initialization(N,dim,ub,lb);

    % L-SHADE options: [population_size, max_FE, 0]
 
    %max_fe = 2000*dim;
    options    = [N, max_fes, 1];

    % Run L-SHADE N_runs times
    fopt_all   = zeros(1, N_runs);
    fopt_all_1 = zeros(1, N_runs);
    fopt_all_2 = zeros(1, N_runs);
    fopt_all_3 = zeros(1, N_runs);
    fopt_all_4 = zeros(1, N_runs);
    fopt_all_5 = zeros(1, N_runs);
    for run = 1 : N_runs
        rng(run + fi*100, 'twister');
        fprintf("run: %s %d times\n", strcat(['F',num2str(fi)]),run);
        
        [Xbest,fbest,curve]  = DE_rand_1(x0,T,N,lb,ub,dim,fobj);
        fopt_all(run)        = fbest;
        
        tic;
        [~, fbest]           = SHADE_call(fobj, lb, ub, dim, max_fes, N, x0);
        times = toc;
        fopt_all_1(1,run)    = fbest; 
        
        tic;
        [Xbest, fbest]       = LSHADE(x0,fobj, lb, ub, dim, max_fes, N);
        times = toc;
        fopt_all_2(1,run)    = fbest;

        tic;
        [~, fbest]           = jSO(fobj, lb, ub, dim, max_fes, N,x0);
        times=toc;
        fopt_all_3(1,run)    = fbest;
      
        
        tic;
        [Xbest,fbest,curve]  = EhdDE(x0, T_ehdde, N, lb, ub, dim, fobj);
        times = toc;
        fopt_all_4(1,run)    = fbest; 
        
        
        tic;
        [fbest,Xbest,curve]  = HODE_para_NoEES(fobj, N, T_hode, lb, ub, dim, x0);
        times = toc;
        fopt_all_5(1,run)    = fbest; 
        
    end
    
    % Statistics
    f_mean(fi)    = mean(fopt_all);
    f_std(fi)     = std(fopt_all);
    f_best(fi)    = min(fopt_all);
    f_worst(fi)   = max(fopt_all);
    
    f_mean_1(fi)  = mean(fopt_all_1);
    f_std_1(fi)   = std(fopt_all_1);
    f_best_1(fi)  = min(fopt_all_1);
    f_worst_1(fi) = max(fopt_all_1);
    
    f_mean_2(fi)  = mean(fopt_all_2);
    f_std_2(fi)   = std(fopt_all_2);
    f_best_2(fi)  = min(fopt_all_2);
    f_worst_2(fi) = max(fopt_all_2);
    
    f_mean_3(fi)  = mean(fopt_all_3);
    f_std_3(fi)   = std(fopt_all_3);
    f_best_3(fi)  = min(fopt_all_3);
    f_worst_3(fi) = max(fopt_all_3);
    
    f_mean_4(fi)  = mean(fopt_all_4);
    f_std_4(fi)   = std(fopt_all_4);
    f_best_4(fi)  = min(fopt_all_4);
    f_worst_4(fi) = max(fopt_all_4);

    f_mean_5(fi)  = mean(fopt_all_5);
    f_std_5(fi)   = std(fopt_all_5);
    f_best_5(fi)  = min(fopt_all_5);
    f_worst_5(fi) = max(fopt_all_5);

%     fprintf('F%-2d  dim=%-3d  mean=%-12.4e  std=%-10.4e  best=%-12.4e\n', ...
%         fi, dim, f_mean(fi), f_std(fi), f_best(fi));
end
elapsed = toc;
fprintf('========================================================\n');
fprintf('Total time: %.2f s\n', elapsed);
fprintf('\n--- Final Results Table ---\n');
fprintf('Func\tMean\t\tStd\t\tBest\n');
data = [f_mean',f_mean_1',f_mean_2',f_mean_3',f_mean_4',f_mean_5',...
        f_std', f_std_1', f_std_2', f_std_3', f_std_4', f_std_5'];
xlswrite(strcat(['cec2022_test_data_rand_run',num2str(max_fes),'_',num2str(N),'_',num2str(D),'.xlsx']),data);
% for fi = func_list
%     fprintf('F%-2d\t%.4e\t%.4e\t%.4e\n', fi, f_mean(fi), f_std(fi), f_best(fi));
% end
% 
% for fi = func_list
%     fprintf('F%-2d\t%.4e\t%.4e\t%.4e\n', fi, f_mean_1(fi), f_std_1(fi), f_best_1(fi));
% end
