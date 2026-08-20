%% 对比不同的DD算法
%% DE start
function [Xbest,fbest,curve] = EhdDE(X,Max_iter,N,lb,ub,dim,fobj)

F     = 2;
fit   = zeros(1,N);
fbest = inf;
curve = zeros(1,Max_iter);
pop   = N;

% factors = Get_Factor_List(si,Max_iter,N,dim);
Xbest = X(1,:);
for t = 1: Max_iter
    
    for i = 1:N
        fit(i) = fobj(X(i,:));
        if fit(i) < fbest
            fbest = fit(i);
            Xbest = X(i,:);
        end
    end
    
    %r2  = factors(t);
    %     if si > 3
    %        r2  = factors(t);
    %     end
    
    r2 = 0.6 * ( exp( (Max_iter/(Max_iter + t)) ) - 1 );
    for i = 1: N
        [~,idx1] = sort(fit);
        supoptimal_X    = X(idx1(2:4),:);
        Index_arr = [1:pop]; %种群下标
        Index_arr(idx1(1:4)) = []; %%去掉适应度前4优种群
        rnd_idx = randperm( pop-4,3); %3个随机下标
        rdi= Index_arr(rnd_idx);%3个随机种群的下标
        %supoptimal_X    = X(idx1(2:4),:);
        
        %       %% 取次优的3个种群，差分，受适应度值次优的前3个种群引导，但只随机取一个维度，并考虑次优秀与全局最优Xbest的信息
        
        %% 次优的3个种群对随机种群rdi(2)和rdi(3)尽心引导,但为了增加种群多样性，只随机取一个维度的值。
        X_top   = ((supoptimal_X(1,:) + supoptimal_X(2,:) + supoptimal_X(3,:))/3 + X(rdi(1),:));
        %rnd_dim = randperm(dim,1);
        %V_rdm   = X_top(1,rnd_dim);
        V_rdm = X_top./2;
        V(1,:)  = V_rdm +    r2 * (X(rdi(3),:) - X(rdi(2),:) ) +  r2 .*( Xbest - supoptimal_X(1,:) );
        
        X_top   = ((supoptimal_X(1,:) + supoptimal_X(2,:) + supoptimal_X(3,:))/3 + X(rdi(2),:));
        %rnd_dim = randperm(dim,1);
        %V_rdm   = X_top(1,rnd_dim);
        V_rdm = X_top./2;
        V1(1,:) = V_rdm  +    r2 .* (X(rdi(3),:) - X(rdi(1),:)) + r2.*( Xbest - supoptimal_X(2,:));
        
        
        X_top   = ((supoptimal_X(1,:) + supoptimal_X(2,:) + supoptimal_X(3,:))/3 + X(rdi(3),:));
        %rnd_dim = randperm(dim,1);
        %V_rdm   = X_top(1,rnd_dim);
        V_rdm = X_top./2;
        V2(1,:) = V_rdm  +   Levy_2(dim,0.5) .* (X(rdi(2),:) - X(rdi(1),:)) + Levy_2(dim,0.5).*(Xbest - supoptimal_X(3,:));
        
        
        rdi = randi([1,pop],1,3);
        X2_new    = X(rdi(1),:) + r2 * (X(rdi(2),:) - X(rdi(3),:));  %%保留ED本身的优势
        Vs = [V;V1;V2;X2_new];
        Vts = [X(i,:);X(i,:);X(i,:);X(i,:)];
        
        
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
        
        
        f2 = fobj(Vts(1,:));
        f3 = fobj(Vts(2,:));
        f4 = fobj(Vts(3,:));
        f6 = fobj(Vts(4,:));
        
        [f5,f_id] = min([f2 f3 f4 f6]);
        if f5 < fit(i)
            fit(i)  = f5;
            X(i,:)  = Vts(f_id,:);
            if f5 < fbest
                Xbest  = Vts(f_id,:);
                fbest  = f5;
            end
        end
        
    end
    %if t > 0.8 * Max_iter
        
        %beta = (2*(1 - t/(Max_iter + eps))^2)*cos(t);
        beta = randn(1,1);
        lv   = Levy(dim);
        %F = 2*(1/t);
        F = 1/t;
        for i = 1 : N
            %X2_new      = X(i,:)+ (1-2*rand(1,1)) .* ( lb.*F + rand(1,1) .* (ub.*F-lb.*F));
            X2_new      = X(i,:)+  (beta).*(lb.*F + rand(1,1) .* (ub.*F-lb.*F));
            X2_new(1,:) = max(X2_new(1,:),lb.*F);
            X2_new(1,:) = min(X2_new(1,:),ub.*F); %% 确定在缩减的区域
            X2_new(1,:) = max(X2_new(1,:),lb);
            X2_new(1,:) = min(X2_new(1,:),ub);
            F_P2        = fobj(X2_new(1,:));
            if F_P2 < fit(i)
                X(i,:) = X2_new(1,:);
                fit(i) = F_P2;
                if F_P2 < fbest
                    fbest  = F_P2;
                    Xbest  = X2_new;
                end
            end
        end
    %end
    
    if t > 0.8 * Max_iter
        %%邻域
        n = 0;
        neibor = zeros(pop,dim);
        if fbest == inf
            numR = round(gen_random(0.1,0.2,1) * pop);
            n = numR;
            %%计算欧式距离
            for i = 1:pop
                dist(i) = sqrt(sum(power((X(i,:)- Xbest(1,:)),2)));
            end
            [~,idx] = sort(dist);
            neibor  = X(idx(1:numR),:);
            %xlswrite("origin_n.xls",neibor);
            %avg = mean(neibor);
            for i = 1:n
                rrr = rand();
                if rrr >= 0.5
                    neibor(i,:) = Xbest + r2*(Xbest - neibor(i,:));%benchmark函数最多
                else
                    neibor(i,:) = Xbest * eps .* Levy(dim); %%附近邻域中最优的X(idx(1),:)
                end
            end
        end
        if fbest ~= inf
            %fprintf("fff\n");
            n = 0;
            for j = 1:pop  %%===通过适应度值构建邻居
                r = gen_random(1,1.3,1);
                %r = K;
                if fit(j) < r * fbest && fit(j) ~= fbest
                    rrr = rand();
                    n = n+1;
                    if rrr > 0.5
                        neibor(i,:) = Xbest + r2*(Xbest - X(j,:));%benchmark函数最多
                    else
                        neibor(n,:) = Xbest  .* Levy(dim) * eps;
                    end
                    
                end
            end
%             if n < (1/10) * pop %% neibor为空
%                 %if n == 0
%                 %找离GBestX最近的50%到70%的，距离最近的点。
%                 numR = round(gen_random(0.1,0.2,1) * pop);
%                 n = numR;
%                 %%计算欧式距离
%                 for i = 1:pop
%                     dist(i) = sqrt(sum(power((X(i,:)- Xbest(1,:)),2)));
%                 end
%                 [~,idx] = sort(dist);
%                 neibor  = X(idx(1:numR),:);
%                 %xlswrite("origin_n.xls",neibor);
%                 %avg = mean(neibor);
%                 for i = 1:n
%                     rrr = rand();
%                     if rrr >= 0.5
%                         neibor(i,:) = Xbest + r2*(Xbest - neibor(i,:));%benchmark函数最多
%                     else
%                         neibor(i,:) = Xbest * eps .* Levy(dim); %%附近邻域中最优的X(idx(1),:)
%                     end
%                     
%                 end
%                 %xlswrite("neibor.xls",neibor);
%             end
        end
        
        FU = neibor(1:n,:) > ub;
        FL = neibor(1:n,:) < lb;
        neibor(1:n,:) = (neibor(1:n,:) .* (~(FU+FL))) + ub .* FU + lb .* FL;
        
        %% 邻域内找最优值和最优种群
        for j = 1: n
            %
            f1 = fobj(neibor(j,:));
            if  f1 < fbest
%                 fprintf("fitness==%.16f,n===%d,t====%d,GbestF=====%.15f\n",f1,n,t,fbest);
                fbest      = f1;
                Xbest      = neibor(j,:);
            end
        end
    end
    curve(t) = fbest;
end