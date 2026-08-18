
fhd = str2func('cec17_func');

%only for D=2,10,30,50,100.

runTimes         = 1;
T                = 2000; 
N                = 100;

lb               = -100;
ub               = 100;
options.max_iter = T;
popsNum          = N;
dim              = 100;

fun_num          = 2;
% 预分配存储变量（避免 parfor 动态扩展报错）
avg_arr = zeros(fun_num, 12);
std_arr = zeros(fun_num, 12);
min_arr = zeros(fun_num, 12);
max_arr = zeros(fun_num, 12);


best_pops_for_problem = cell(fun_num, 1);

popinfo = struct('lb', lb, ...
    'ub', ub, ...
    'dim', dim);
x0 = initialization( N, dim, ub, lb);

All_data = cell(fun_num,1);

for j = 1: 6
    
    curves = zeros(12, T);
    fprintf('function is F%d\n',j);
    
%     fun_name = strcat(['F',num2str(j)]);
%     [lb,ub,dim,fobj]    = Get_Functions_details(strcat(fun_name));
     

    fun_name = j;
    [fobj,dim,lb,ub]    = get_function(fun_name);    
    
     popInfo = struct('lb', lb, ...
        'ub',  ub, ...
        'dim', dim);
    x0 = initialization(N,dim,ub,lb);
    
    lowerbound = lb;
    upperbound = ub;
    dimension  = dim;
    pops       = x0;
    % 预分配存储（必须 inside parfor，每个 worker 独立）
    gbest_arr_1  = zeros(1, runTimes);
    gbest_arr_2  = zeros(1, runTimes);
    gbest_arr_3  = zeros(1, runTimes);
    gbest_arr_4  = zeros(1, runTimes);
    gbest_arr_5  = zeros(1, runTimes);
    gbest_arr_6  = zeros(1, runTimes);
    gbest_arr_7  = zeros(1, runTimes);
    gbest_arr_8  = zeros(1, runTimes);
    gbest_arr_9  = zeros(1, runTimes);
    gbest_arr_10 = zeros(1, runTimes);
    gbest_arr_11 = zeros(1, runTimes);
    gbest_arr_12 = zeros(1, runTimes);
    gbest_arr_13 = zeros(1, runTimes);
    
    best_param_1 = NaN(runTimes, dim);
    best_param_2 = NaN(runTimes, dim);
    best_param_3 = NaN(runTimes, dim);
    best_param_4 = NaN(runTimes, dim);
    best_param_5 = NaN(runTimes, dim);
    best_param_6 = NaN(runTimes, dim);
    best_param_7 = NaN(runTimes, dim);
    best_param_8 = NaN(runTimes, dim);
    best_param_9 = NaN(runTimes, dim);
    best_param_10 = NaN(runTimes, dim);
    best_param_11 = NaN(runTimes, dim);
    best_param_12 = NaN(runTimes, dim);
    
    % --- 多次运行 ---
    for i = 1: runTimes
        fprintf('\t F%d runtimes is %d\n', j, i);
        
        % 每次设置不同随机种子，避免并行重复
        rng(i + j*100, 'twister');
        
        % === 调用你的算法们 ===
        tic;
        %         %fhd,SearchAgents,Max_iterations,lowerbound,upperbound,dimension,pops
        [fbest,Xbest,curve]     = HODE_para_NoEES(fobj,N,T,lb,ub,dim,x0);
        times               = toc;
        runtims_arr_1(j)    = times;
        gbest_arr_1(i)      = fbest;
        best_param_1(i,:)   = Xbest;
        curves(1,:)         = curve; 
        
        tic;
        [Xbest,fbest,curve] = EnhancedDE_TSP_BKNP(x0,T,N,lb,ub,dim,fobj,3);
        times = toc;
        runtims_arr_2(j)    = times;
        gbest_arr_2(i)      = fbest;
        best_param_2(i,:)   = Xbest;
        curves(2,:)         = curve;
         
        %
        tic;
        [Xfood, fval,~,~,~, ~,curve] = ISO(N,T,lb,ub,dim,fobj,x0);
        times = toc;
        runtims_arr_3(j)    = times;
        gbest_arr_3(i)      = fval;
        best_param_3(i,:)   = Xfood;
        curves(3,:)         = curve;
        %
        tic;
        [Xbest,fbest,curve] = DE_standard(x0,T,N,lb,ub,dim,fobj);
        times = toc;
        runtims_arr_4(j)    = times;
        gbest_arr_4(i)      = fbest;
        best_param_4(i,:)   = Xbest;
        curves(4,:)         = curve;
         %
        % %
        tic;
        [Best_score,Best_pos,curve]= OOBO(fobj,N,T,lb,ub,dim,x0);
        times = toc;
        runtims_arr_5(j)    = times;
        gbest_arr_5(i)      = Best_score;
        best_param_5(i,:)   = Best_pos;
        curves(5,:)         = curve;
        
        tic;
        [~,Xbest, curve]        = IVY(N,T,lb,ub,dim,fobj,x0);
        times = toc;
        runtims_arr_6(j)    = times;
        gbest = min(curve);
        gbest_arr_6(i)      = gbest;
        best_param_6(i,:)   = Xbest;
        curves(6,:)         = curve;
        
        
        tic;
        [Best_score,Xbest,curve] = WOA(fobj,N,T,lb,ub,dim,x0);
        times = toc;
        runtims_arr_7(j)    = times;
        gbest_arr_7(i)      = Best_score;
        best_param_7(i,:)   = Xbest;
        curves(7,:)         = curve;
       
        
%         tic;
%         [Best_score,Best_SW,curve] = SWO( N,T,ub,lb,dim,fobj,x0,popInfo);
%         times = toc;
%         runtims_arr_8(j)    = times;
%         gbest_arr_8(i)      = Best_score;
%         best_param_8(i,:)   = Best_SW;
%         curves(8,:)         = curve;
%         
        
        tic
        [Best_score,Xbest,curve]      = HO(N,T,lb,ub,dim,fobj,x0,popInfo);
        times = toc;
        runtims_arr_8(j)    = times;
        gbest_arr_8(i)      = Best_score;
        best_param_8(i,:)   = Xbest;
        curves(8,:)         = curve;
    
        %
        tic;
        [Xbest, Food_Score,curve]  = HBA(fobj, dim,lb,ub,T,N,x0,popInfo);
        times = toc;
        runtims_arr_9(j)    = times;
        gbest_arr_9(i)      = Food_Score;
        best_param_9(i,:)   = Xbest;
        curves(9,:)         = curve;

        %
        %
        tic;
        [Xbest,fvalbest, curve] = BWO(N,T,lb,ub,dim,fobj,x0);
        times = toc;
        runtims_arr_10(j)    = times;
        gbest_arr_10(i)      = fvalbest;
        best_param_10(i,:)   = Xbest;
        curves(10,:)         = curve;

        %
        tic;
        [Best_score,Xbest,curve] = GOA( N,T,lb,ub,dim,fobj,x0);
        times = toc;
        runtims_arr_11(j)    = times;
        gbest_arr_11(i)      = Best_score;
        best_param_11(i,:)   = Xbest;
        curves(11,:)         = curve;

        %Best_score
        
        %
        tic;
        [Best_F,Xbest,curve]       = RSA( N,T,lb,ub,dim,fobj,x0);
        times = toc;
        runtims_arr_12(j)    = times;
        gbest_arr_12(i)      = Best_F;
        curves(12,:)         = curve;
        best_param_12(i,:)   = Xbest;
       
        
    end
    %xlswrite(strcat(['eng/curve/curve_F',num2str(j),'_',num2str(T),'_',num2str(N),'.xlsx']),curves);
    
    All_data{j,:}=[gbest_arr_1',gbest_arr_2',gbest_arr_3',gbest_arr_4',gbest_arr_5',...
                   gbest_arr_6',gbest_arr_7',gbest_arr_8',gbest_arr_9',gbest_arr_10',...
                   gbest_arr_11',gbest_arr_12'];   
    
    % --- 保存每个函数的最佳个体 ---
        all_best_pop = [best_param_1; best_param_2; best_param_3; best_param_4; ...
            best_param_5; best_param_6; best_param_7; best_param_8; ...
            best_param_9; best_param_10; best_param_11;best_param_12];
        best_pops_for_problem{j,1} = all_best_pop;
    
    
    % --- 统计结果 ---
    avg_arr(j,:) = [mean(gbest_arr_1), mean(gbest_arr_2), mean(gbest_arr_3), mean(gbest_arr_4), ...
        mean(gbest_arr_5), mean(gbest_arr_6), mean(gbest_arr_7), mean(gbest_arr_8), ...
        mean(gbest_arr_9), mean(gbest_arr_10), mean(gbest_arr_11),mean(gbest_arr_12)];
    std_arr(j,:) = [std(gbest_arr_1), std(gbest_arr_2), std(gbest_arr_3), std(gbest_arr_4), ...
        std(gbest_arr_5), std(gbest_arr_6), std(gbest_arr_7), std(gbest_arr_8), ...
        std(gbest_arr_9), std(gbest_arr_10), std(gbest_arr_11),std(gbest_arr_12)];
    min_arr(j,:) = [min(gbest_arr_1), min(gbest_arr_2), min(gbest_arr_3), min(gbest_arr_4), ...
        min(gbest_arr_5), min(gbest_arr_6), min(gbest_arr_7), min(gbest_arr_8), ...
        min(gbest_arr_9), min(gbest_arr_10), min(gbest_arr_11),min(gbest_arr_12)];
    max_arr(j,:) = [max(gbest_arr_1), max(gbest_arr_2), max(gbest_arr_3), max(gbest_arr_4), ...
        max(gbest_arr_5), max(gbest_arr_6), max(gbest_arr_7), max(gbest_arr_8), ...
        max(gbest_arr_9), max(gbest_arr_10), max(gbest_arr_11),max(gbest_arr_12)];
    
end
% runTimes = [runtims_arr_1',runtims_arr_2',runtims_arr_3',runtims_arr_4',runtims_arr_5',...
%             runtims_arr_6',runtims_arr_7',runtims_arr_8',runtims_arr_9',runtims_arr_10',...
%             runtims_arr_11',runtims_arr_12'];
% xlswrite(strcat(['cec_23_bnm/cec_bnch_runtimes.xls']), runTimes);

%% for boxfig
% for i = 1: fun_num
%     data1 = All_data{i,:};
%     xlswrite(strcat(['cec_eng/boxfig/eng_funs_for_boxfig_F',num2str(i),'_',num2str(T),'_',num2str(N),'.xlsx']), data1);
% end
% % for statistic data 

data = [avg_arr, std_arr, min_arr, max_arr];
xlswrite(strcat(['cec_eng/best_test_data_',num2str(T),'_',num2str(N),'.xlsx']), data);

%% for best param / need set runtimes = 1
for i = 1: fun_num
    x = best_pops_for_problem(i,1);
    xlswrite(strcat(['cec_eng/best_param_F',num2str(i),'.xlsx']), cell2mat(x));
end


