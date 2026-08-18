
%%% Designed and Developed by Mohammad Dehghani %%%

function[Best_score,Best_pos,OOBO_curve]=HODE_para_NoEES(fobj,SearchAgents,Max_iterations,lowerbound,upperbound,dimension,pops)

lowerbound=ones(1,dimension).*(lowerbound);                              % Lower limit for variables
upperbound=ones(1,dimension).*(upperbound);                              % Upper limit for variables

dim = dimension;
lb = lowerbound;
ub = upperbound;
%%
% for i=1:dimension
%     X(:,i) = lowerbound(i)+rand(SearchAgents,1).*(upperbound(i) - lowerbound(i));                          % Initial population
% end
X1 = pops;
X2 = pops;
%F = 2;
pop = SearchAgents;
fit1 = zeros(1,SearchAgents);
fit2 = zeros(1,SearchAgents);
for i =1:SearchAgents
    L1 = X1(i,:);
    L2 = X2(i,:);
    %fitness=feval(fhd,Positions(i,:)',varargin{:});
    fit1(i) = fobj(L1(1,:));
    fit2(i) = fobj(L2(1,:));
end
%%
fbest2 = inf;
fbest1 = inf;
fbest  = inf;
Max_iter = Max_iterations;

Xbest = X1(1,:);
Xbest1 = Xbest;
Xbest2 = Xbest;
for t=1:Max_iterations
    
    %Curve_2dim(t,:) = [X(1,1),X(2,1)];
    %     if t == 10
    %         xlswrite(strcat(['Bnchmark_pop_dis_',func2str(fobj),'.xlsx']),X(:,1:2));%%存2个维度就够了
    %     end
    %% update the global best (fbest)
    [best1 , location] = min(fit1);
    if t==1
        Xbest1 = X1(location,:);                                           % Optimal location
        fbest1 = best1;                                           % The optimization objective function
    elseif best1 < fbest1
        fbest1 = best1;
        Xbest1 = X1(location,:);
    end
    
    %r2 = Get_Factor_List( sch,Max_iter,N,dim)
    %% update agents position
    K=1:SearchAgents; %Eq(9)
    
    %% OK
    
    %r2 = 0 - (0 - 2) * rand() + 0.5*randn(1,1); %% 11
    for i=1:SearchAgents
        %         rr = rand();
        %         if rr <=0.5
        k = randperm(SearchAgents-i+1,1);
        k = K(k);
        
        if k==i
            k=k+1;
            if k>SearchAgents
                k=k-2;
            end
        end
        %tt = find(K==k)
        K( find(K==k) ) = [];
        X_k = X1(k,:);
        I=round(1+rand); %Eq(12)
        
        if fit1(i) > fit1(k)
            X_new1 = X1(i,:)+ rand(1,dimension).*(X_k - I .*  X1(i,:)); %Eq(11)
        else
            X_new1 = X1(i,:)+ rand(1,dimension).*(X1(i,:) - X_k); %Eq(11)
        end
        X_new1 = max(X_new1,lowerbound);
        X_new1 = min(X_new1,upperbound);
        
        % based on Eq(13)
        ff1 = fobj(X_new1(1,:));
        if ff1 <= fit1 (i)
            X1(i,:) = X_new1;
            fit1 (i)= ff1;
        end
        
        
        %% DE start
        %%r2  = 2- 2/(1 + exp(-t)) + eps;
        %r2 = 2 * cos(((1 - 30 * t / Max_iter))^5);
        %r2 = 2 * (1-t/Max_iter)^2 *cos(t);
        %r2 = 2 * (1-t/Max_iter)^(2*(t/Max_iter));
    end
    
    r2 = 0.6 * ( exp( (Max_iter/(Max_iter + t)) ) - 1 ); %%26
%     [best2 , location2] = min(fit2);
%     if t==1
%         Xbest2 = X2(location2,:);                                           % Optimal location
%         fbest2 = best2;                                           % The optimization objective function
%     elseif best2 < fbest2
%         fbest2 = best2;
%         Xbest2 = X2(location2,:);
%     end
    for i = 1:pop
        fit2(i) =  fobj(X2(i,:));
        
        if fit2(i) < fbest2
            fbest2 = fit2(i);
            Xbest2 = X2(i,:);
        end
    end
    
    for i = 1:pop
        [~,idx1] = sort(fit2);
        supoptimal_X    = X2(idx1(2:4),:);
        Index_arr = [1:pop]; %种群下标
        Index_arr(idx1(1:4)) = []; %%去掉适应度前4优种群
        rnd_idx = randperm( pop-4,3); %3个随机下标
        rdi= Index_arr(rnd_idx);%3个随机种群的下标
        %supoptimal_X    = X(idx1(2:4),:);
        
        %       %% 取次优的3个种群，差分，受适应度值次优的前3个种群引导，但只随机取一个维度，并考虑次优秀与全局最优Xbest的信息
        
        %% 次优的3个种群对随机种群rdi(2)和rdi(3)尽心引导,但为了增加种群多样性，只随机取一个维度的值。
        X_top   = ((supoptimal_X(1,:) + supoptimal_X(2,:) + supoptimal_X(3,:))/3 + X2(rdi(1),:));
        rnd_dim = randperm(dim,1);
        %V_rdm   = X_top(1,rnd_dim);
        V_rdm   = X_top./2;
        V(1,:)  = V_rdm +    r2 * (X2(rdi(3),:) - X2(rdi(2),:) ) +  r2 .*( Xbest2 - supoptimal_X(1,:) );
        
        X_top   = ((supoptimal_X(1,:) + supoptimal_X(2,:) + supoptimal_X(3,:))/3 + X2(rdi(2),:));
        rnd_dim = randperm(dim,1);
        %V_rdm   = X_top(1,rnd_dim);
        V_rdm = X_top./2; 
        V1(1,:) = V_rdm  +   r2 * (X2(rdi(3),:) - X2(rdi(1),:)) + r2 .*( Xbest2 - supoptimal_X(2,:));
        
        
        X_top   = ((supoptimal_X(1,:) + supoptimal_X(2,:) + supoptimal_X(3,:))/3 + X2(rdi(3),:));
        rnd_dim = randperm(dim,1);
        %V_rdm   = X_top(1,rnd_dim);
        V_rdm = X_top./2;
        V2(1,:) = V_rdm  +   Levy_2(dim,0.5) .* (X2(rdi(2),:) - X2(rdi(1),:)) +  Levy_2(dim,0.5) .*(Xbest2 - supoptimal_X(3,:));
        
        %Vs = [V;V1;V2];
        
        rdi = randi([1,pop],1,3);
        X2_new    = X2(rdi(1),:) + r2 * (X2(rdi(2),:) - X2(rdi(3),:));  %%保留ED本身的优势
        Vs = [V;V1;V2;X2_new];
        
        Vts = [X2(i,:);X2(i,:);X2(i,:);X2(i,:)];
        
        for k = 1:4
            for d = 1 : dim
                r3 = rand();
                dimr = randi([1,dim],1);
                if r3 < 0.9 || dimr == d
                    Vts(k,d) = Vs(k,d);
                end
            end
        end
        
        FU = Vts > ub;
        FL = Vts < lb;
        Vts = Vts .* (~(FU+FL)) + ub .* FU + lb .* FL;
        
        %% 3个中去最小 好于3个均值
        f2 = fobj(Vts(1,:));
        f3 = fobj(Vts(2,:));
        f4 = fobj(Vts(3,:));
        f6 = fobj(Vts(4,:));
       
        
        [f5,f_id] = min([f2 f3 f4 f6]);
        if f5 < fit2(i)
            fit2(i)  = f5;
            X2(i,:)  = Vts(f_id,:);
            if f5 < fbest2
                Xbest2  = Vts(f_id,:);
                fbest2  = f5;
            end
        end
    end
    %% Local search,每个种群都要局部搜索
    
    beta = randn(1,1);
    for i=1 : SearchAgents
        %X2_new      = X2(i,:)+ (1-2*rand(1,1)) .* ( lb./log(1+t) + rand(1,1).*(ub./log(1+t)-lb./log(1+t)));
        X2_new      = X2(i,:)+ (beta) .* ( lb./t + rand(1,1) .* (ub./t-lb./t));
        %X2_new      = X2(i,:)+ ( Levy(dim)) .* ( lb./log(1+t) + Levy(dim).*(ub./log(1+t)-lb./log(1+t)));
        X2_new(1,:) = max(X2_new(1,:),lb./t);
        X2_new(1,:) = min(X2_new(1,:),ub./t);
        X2_new(1,:) = max(X2_new(1,:),lb);
        X2_new(1,:) = min(X2_new(1,:),ub); 
        
        % update position based on Eq (7)
        %L = X2_new(i,:);
        F_P2 = fobj(X2_new(1,:));
        if F_P2 < fit2(i)
            X2(i,:) = X2_new(1,:);
            fit2(i) = F_P2;
            if F_P2 < fbest2
                fbest2  = F_P2;
                Xbest2  = X2_new;
            end
        end
    end % END for i=1:SearchAgents
    %
    %% neighborhood  start
    if t > Max_iterations*0.8
        n = 0; 
        neibor = zeros(pop,dim);
        if fbest2 == inf
            numR = round(gen_random(0.1,0.2,1) * pop);
            n = numR;
            %%计算欧式距离
            for i = 1:pop
                dist(i) = sqrt(sum(power((X2(i,:)- Xbest2(1,:)),2)));
            end
            [~,idx] = sort(dist);
            neibor  = X2(idx(1:numR),:);
            %xlswrite("origin_n.xls",neibor);
            %avg = mean(neibor);
            for i = 1:n
                rrr = rand();
                if rrr >= 0.5
                    %neibor(i,:) = GBestX + rand() * abs(GBestX - neibor(i,:)).* Levy(dim);%benchmark函数最多
                    neibor(i,:) = Xbest2 +  r2*(Xbest2 - neibor(i,:));
                else
                    neibor(i,:) = Xbest2 * eps .* Levy(dim); %%附近邻域中最优的X(idx(1),:)
                end
            end
        end
        if fbest2 ~= inf
            %fprintf("fff\n");
            n = 0;
            for j = 1:pop  %%===通过适应度值构建邻居
                r = gen_random(1,1.3,1);
                %r = K;
                if (fit2(j) < r * fbest2) && (fit2(j) ~= fbest2)
                    rrr = rand();
                    n = n+1;
                    if rrr > 0.5
                        %neibor_Pops(i,:) = pops(idx(1,i),:) + T_distribution(2,0.5,popInfo.dim) .* eps
                        neibor(n,:) = Xbest2  + r2 * (Xbest2 - X2(j,:)); %benchmark函数最多
                        %                     r2 = F * (1- t/Max_iter)^(2);
                        %                     neibor(i,:) = GBestX + r2*cos(t)*sin(t)*abs( median(neibor) - neibor(i,:));


                    else
                        neibor(n,:) = Xbest2  .* Levy(dim) * eps;
                    end

                end
            end
            if n < (1/10) * pop %% neibor为空
                %if n == 0
                %找离GBestX最近的50%到70%的，距离最近的点。
                numR = round(gen_random(0.1,0.2,1) * pop);
                n = numR;
                %%计算欧式距离
                for i = 1:pop
                    dist(i) = sqrt(sum(power((X2(i,:)- Xbest2(1,:)),2)));
                end
                [~,idx] = sort(dist);
                neibor  = X2(idx(1:numR),:);
                %xlswrite("origin_n.xls",neibor);
                %avg = mean(neibor);
                for i = 1:n
                    rrr = rand();
                    if rrr >= 0.5
                        neibor(i,:) = Xbest2 + r2*(Xbest2 - neibor(i,:));%benchmark函数最多
                        %                    r2 = F * (1- t/Max_iter)^(2);
                        %                    neibor(i,:) = GBestX + r2*cos(t)*sin(t)*abs(median(neibor) - neibor(i,:));

                    else
                        neibor(i,:) = Xbest2 * eps .* Levy(dim); %%附近邻域中最优的X(idx(1),:)
                    end

                end
                %xlswrite("neibor.xls",neibor);
            end
        end

        FU = neibor(1:n,:) > ub;
        FL = neibor(1:n,:) < lb;
        neibor(1:n,:) = (neibor(1:n,:) .* (~(FU+FL))) + ub .* FU + lb .* FL;

        %% 邻域内找最优值和最优种群
        for j = 1: n
            %
            f1 = fobj(neibor(j,:));
            if  f1 < fbest2
%                 fprintf("HODE_no_CH_fitness==%.16f,n===%d,t====%d,GbestF=====%.15f\n",f1,n,t,fbest2);
                fbest2      = f1;
                Xbest2      = neibor(j,:);
            end
        end
    end
    
    if fbest1 < fbest2
        fbest = fbest1;
        Xbest = Xbest1;
    else
        fbest = fbest2;
        Xbest = Xbest2;
    end
    %if fbest2 < fbest
       
    %end
    OOBO_curve(t) = fbest;
end
%xlswrite(strcat(['data\bnch_2dim_iteration_',func2str(fobj),'.xls']),Curve_2dim);
Best_score = fbest;
Best_pos   = Xbest;
%OOBO_curve=best_so_far;
end

