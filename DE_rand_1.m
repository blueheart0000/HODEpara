%% 对比不同的DD算法
%% DE start
%% DE rand/1 == 标准DE


function [Xbest,fbest,curve] = DE_rand_1(X,Max_iter,N,lb,ub,dim,fobj)

fit = zeros(N,size(X,2));
fbest = inf;
F = 2;
CR = 0.9;
curve = zeros(1,Max_iter);

pop = N;

for t = 1: Max_iter
    for i = 1:N
        fit(i) = fobj(X(i,:));
        if fit(i) < fbest
            fbest = fit(i);
            Xbest = X(i,:);
        end
    end
    for i = 1: N
        %r2  = F * (1- (t/Max_iter))^2 * cos(t);
        
        %r2 = 2 * cos(((1 - 30 * t / Max_iter))^5);
        rdi = randperm(pop,3);
        F = rand();
        V = X(rdi(1),:) + F *( X(rdi(2),:) - X(rdi(3),:));
        X_new = X(i,:);
        for d = 1 : dim
            dimr = randi([1,dim],1);
            r3 = rand();
            if r3 < 0.9 || dimr == d
                X_new(1,d) = V(1,d);
            end
        end
        
        
        FU = X_new > ub;
        FL = X_new < lb;
        X_new = X_new .* (~(FU+FL)) + ub .* FU + lb .* FL;
        
        %% 3个中去最小 好于3个均值
        f2 = fobj(X_new(1,:));
        
        if f2 < fit(i)
            fit(i)  = f2;
            X(i,:)  = X_new(1,:);
            if f2 < fbest
                Xbest  = X_new;
                fbest  = f2;
            end
        end
        curve(t) = fbest;
    end
end