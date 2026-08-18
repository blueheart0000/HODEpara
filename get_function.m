function [fobj,dim,lb,ub]= get_function(FunIndex) %��ȡ��Χ

switch FunIndex
    case 1
        fobj = @welded_beam_obj;
        dim = 4;
        lb = [0.125, 0.1, 0.1, 0.1];
        ub = [2, 10, 10, 2];
        
    case 2
        fobj = @spring_obj;
        dim = 3;
        lb = [0.05, 0.25, 2];
        ub = [2.0,  1.3,  15];
        
    case 3
        fobj = @pressure_vessel_obj;
        dim = 4;
        lb = [0, 0, 10, 10];
        ub = [99, 99, 200, 200];
        
    case 4
        fobj = @truss_obj;
        dim = 2;
        lb = [0,0];
        ub = [1,1];
        
    case 5
        fobj = @Speed_reducer_obj;
        dim = 7;
        lb = [2.6, 0.7, 17, 7.3, 7.3, 2.9, 5];
        ub = [3.6, 0.8, 28, 8.3, 8.3, 3.9, 5.5];
    
    case 6
        fobj = @design_cantilever_beam_obj;
        dim = 5;
        lb = [0.01,0.01,0.01,0.01,0.01];
        ub = [100, 100, 100, 100, 100];
end 
end

%% ����������Լ�� ��һ�� 
function cost = welded_beam_obj(x)   %Լ��
    % ԭʼĿ�꺯����Լ��
    cost = 65535;
    [f, c, ~] = welded_beam_problem(x);
    tol = 0;
    feas = all(c <= tol);
    if feas ~= 0
       cost = f;
    end
%     [raw_cost, constraints] = welded_beam_design(x);
%        if constraints(1)> 0 || constraints(2)> 0   || constraints(3)> 0 || constraints(4)> 0 || constraints(5)> 0  || constraints(6)> 0
%            x = zeros(1,4);
%            cost = inf;
%            %fprintf("cost====%.9f\n",1);
%       
%        end 
%        if  constraints(1) <= 0 && constraints(2) <= 0 && constraints(3) <= 0 && constraints(4)<= 0 && constraints(5)<= 0 && constraints(6)<= 0 
%          
%            cost = raw_cost;
%            %  fprintf("cost====%.9f\n",2);
%        end 
       %fprintf("cost====%.9f\n",cost);

end

%% �����ɣ�Լ�� �ڶ���
function weight = spring_obj(x)  
    % ԭʼĿ�꺯����Լ��
    weight = 65535;
    tol = 0;
    [f, c, ~] = spring_design_problem(x);
    feas = all(c <= tol);
    if feas ~= 0
        weight = f;
    end
    
%     [raw_weight, constraints] = spring_design(x);
%        if  constraints(1)> 0 || constraints(2)> 0 || constraints(3) > 0 || constraints(4) > 0
%            x = zeros(1,3);
%            weight = inf;
%        end 
%        if  constraints(1) <= 0 && constraints(2) <= 0 || constraints(3)<=0 || constraints(4)<=0
%            weight = raw_weight;
%        end
 
end


%% ѹ������ ������ ������
function cost = pressure_vessel_obj(x)
    cost = inf;
    [raw_cost, constraints] = pressure_vessel_design(x);
       if constraints(1)> 0 || constraints(2)> 0   || constraints(3)> 0 || constraints(4)> 0
           x = zeros(1,4);
           cost = 65535;
      
       end 
       if  constraints(1) <= 0 && constraints(2) <= 0 && constraints(3) <= 0 && constraints(4)<= 0
         
           cost = raw_cost;
       end  
    
end

%% �������������� ���ĸ�
function weight = truss_obj(x)
    weight  = inf;
    [raw_weight, constraints] = truss_design(x);
    if constraints(1)> 0 || constraints(2)> 0 || constraints(3)> 0
        weight = 65535;
    end 
    if  constraints(1) <= 0 && constraints(2) <= 0 && constraints(3) <= 0 
        weight = raw_weight;
    end
end

%% ��������� �����
function weight = Speed_reducer_obj(x)
    [raw_weight, constraints] = Speed_reducer(x);
    weight = inf;
    if  constraints(1)> 0 || constraints(2)> 0 || constraints(3)> 0 || constraints(4)> 0 || constraints(5)> 0  || constraints(6)> 0  || constraints(7)> 0  || constraints(8)> 0 ...
            || constraints(9)> 0  || constraints(10)> 0  || constraints(11)> 0
        
           weight = 65535;
    end 
    if  constraints(1) <= 0 && constraints(2) <= 0 && constraints(3) <= 0 && constraints(4)<= 0 && constraints(5)<= 0 && constraints(6)<= 0 && constraints(7)<= 0 ...
            && constraints(8)<= 0 && constraints(9)<= 0 && constraints(10)<= 0 && constraints(11)<= 0
         
           weight = raw_weight;
    end
end


function rs =  design_cantilever_beam_obj(x)
    rs = inf;
    [fitness, constraints] = design_cantilever_beam(x);
    if constraints(1) >0
        rs = 65535;
    else
        rs = fitness;
    end
end





%% ԭʼ���⺯��������֮ǰ�Ķ���һ�£��˴��������·���У�
%===========������������� 
function [cost, constraints] = welded_beam_design(x)
    % ��������
    x1 = x(1); % ����߶� (0.1 <= h <= 2)
    x2 = x(2); % ���쳤�� (0.1 <= l <= 10)
    x3 = x(3); % ���ĺ�� (0.1 <= t <= 10)
    x4 = x(4); % ���Ŀ�� (0.1 <= b <= 2)
    
    % Ŀ�꺯�� (�ܳɱ�)
    cost = 1.10471*x1^2*x2 + 0.04811*x3*x4*(14.0 + x2);
    
    % Լ������
    P = 6000;   % �غ� (lb)
    L = 14;     % ���ĳ��� (in)
    E = 30e6;   % ����ģ�� (psi)
    G = 12e6;   % ����ģ�� (psi)
    tau_max = 13600;   % ����������Ӧ�� (psi)
    sigma_max = 30000; % �����������Ӧ�� (psi)
    delta_max = 0.25;  % �������˲��Ӷ� (in)
    
     % �м�������
    M = P*(L + x2/2);                  % ���
    R = sqrt( (x2^2)/4 + ((x1 + x3)/2)^2 );    % ��������
    J = 2*(sqrt(2*x1*x2)*( x2^2/12 + (x1 + x3)^2/4)); % �����Ծ�  %%
    tau_prime = P/sqrt(2*x1*x2);       % ֱ�Ӽ���Ӧ��  %%�����еĲ�ͬ
    tau_double_prime = (M*R)/J;        % ��������Ӧ��
    tau = sqrt(tau_prime^2 + tau_prime*tau_double_prime*x2/R + tau_double_prime^2); % �ܼ���Ӧ��
    sigma = 6*P*L/(x3^2*x4);            % ����Ӧ��
    P_c   = (4.013*E*sqrt((x3^2*x4^6)/36))/(L^2) * (1 - x3/(2*L)*sqrt(E/(4*G))); % �����غ�
    delta = 6*P*L^3/(E*x3^3*x4);        % �˲��Ӷ�
    
    % ����ʽԼ�� (<=0)
    constraints    = zeros(6,1);
    constraints(1) = tau - tau_max;   % ����Ӧ��Լ��
    constraints(2) = sigma - sigma_max; % ����Ӧ��Լ��
    constraints(3) = x1 - x4;           % ����Լ�� (b >= h)
    %constraints(4) = delta - delta_max; % �Ӷ�Լ��
    constraints(4) = 0.10471*x1^2 + 0.04811*x3*x4*(14.0+x2)-5;
    constraints(5) = delta - delta_max;        % �����غ�Լ��
    constraints(6) = P-P_c;
    %constraints(7) = 1.0471* (x1^2) + 0.04811*x3*x4*(14.0 + x2) - 5.0;
    
end

% ==========�����������
function [weight, constraints] = spring_design(x)
    % ��������
    x1 = x(1); % �߾� (0.05 <= d <= 2.0)
    x2 = x(2); % ��Ȧֱ�� (0.25 <= D <= 1.3)
    x3 = x(3); % ��Ч��Ȧ�� (2 <= N <= 15)
    
    % Ŀ�꺯�� (������С��)
    weight = (x3 + 2)*x2*x1^2;
        
    % ����ʽԼ�� (<=0)
    constraints    = zeros(4,1);
    constraints(1) = 1 - (x3*x2^3)/(71785*x1^4);          
    constraints(2) = (4*(x2^2) - (x1*x2))/(12566*(x2*(x1^3) - x1^4)) + 1/(5108*(x1^2)) - 1; %rsa������ 
    constraints(3) = 1-(140.45*x1)/(x2^2*x3); 
    constraints(4) = (x1 + x2)/1.5 - 1;
end


% ========== ѹ������������� ===========
function [cost, constraints] = pressure_vessel_design(x)
    % ��������
    x1 = x(1); % ������ (0.625 <= Ts <= 1.0)
    x2 = x(2); % ͷ��� (0.625 <= Th <= 1.0)
    x3  = x(3); % �ڰ뾶 (10 <= R <= 200)
    x4  = x(4); % ���� (10 <= L <= 200)
    
    % Ŀ�꺯�����ܳɱ���
    cost = 0.6224*x1*x3*x4 + 1.7781*x2*x3^2 + 3.1661*x1^2*x4 + 19.84*x1^2*x3;
    
    % Լ������
    constraints = zeros(4,1);
    constraints(1) = -x1 + 0.0193*x3;
    constraints(2) = -x2 + 0.00954*x3;
    constraints(3) = -pi*x3^2*x4 - (4/3)*pi*x3^3 + 750*1728;
    constraints(4) = x4 - 240;
   
    
end

% ========== �������������� ===========
function [weight, constraints] = truss_design(x)    %%������
    % ��������
    x1 = x(1); % ��1��� (0.1 <= A1 <= 1)
    x2 = x(2); % ��2��� (0.1 <= A2 <= 1)
    
    % ���ϲ���
    P = 2000;   % �غ� (N)
    l = 100;
    
    % Ŀ�꺯������������
    weight = (2*sqrt(2)*x1 + x2) * l;     %%x1���ģ������в�ͬ
    

    % Լ������
    constraints = zeros(3,1);
    constraints(1) = P*(sqrt(2)*x1 + x2)/(sqrt(2)*x1^2 + 2*x1*x2) - 2000;  % Ӧ��Լ��
    constraints(2) = P*x2/(sqrt(2)*x1^2 + 2*x1*x2)- 2000;
    constraints(3) = P/(sqrt(2)*x2 + x1) - 2000;
    
end

% =========== ������������� =============
function [weight, constraints] = Speed_reducer(x)
    % �������� [b, m, z, l1, l2, d1, d2]
    x1 = x(1);   % �ݿ� (2.6 <= b <= 3.6)
    x2 = x(2);   % ģ�� (0.7 <= m <= 0.8)
    x3 = x(3);   % ���� (17 <= z <= 28)
    x4 = x(4);  % ��1���� (7.3 <= l1 <= 8.3)
    x5 = x(5);  % ��2���� (7.3 <= l2 <= 8.3)
    x6 = x(6);  % ��1ֱ�� (2.9 <= d1 <= 3.9)
    x7 = x(7);  % ��2ֱ�� (5.0 <= d2 <= 5.5)
    
    % Ŀ�꺯������������
    weight = 0.7854*x1*(x2^2)*(3.3333*(x3^2) + 14.9334*x3 - 43.0934) - 1.508*x1*(x6^2 + x7^2) + 7.4777*(x6^3 + x7^3) + 0.7854*(x4*(x6^2) + x5*(x7^2));
    
    % Լ������
    constraints = zeros(11,1);
    
    
    constraints(1) = 27/(x1*x2^2*x3) - 1;
    constraints(2) = 397.5/(x1*x2^2*x3^2) - 1;
    constraints(3) = (1.93*x4^3)/(x2*x3*x6^4) - 1;
    constraints(4) = (1.93*x5^3)/(x2*x3*x7^4) - 1;
    constraints(5) = sqrt((745*x4/(x2*x3))^2 + 16.9e6)/(110*x6^3) - 1;
    constraints(6) = sqrt((745*x5/(x2*x3))^2 + 157.5e6)/(85*x7^3) - 1;
    constraints(7) = (x2*x3)/40 - 1;
    constraints(8) = (5*x2)/x1 - 1;
    constraints(9) = x1/(12*x2) - 1;    %%%�е�����
    
    constraints(10) = ((1.5* x6 )+ 1.9 ) / x4 - 1;
    constraints(11) = ((1.1* x7 )+ 1.9 ) / x5 - 1;
end



function [fitness,constraints] = design_cantilever_beam(x)
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    x5 = x(5);
    
    fitness = 0.0624 * (x1+x2+x3+x4+x5);
    constraints(1) = 61/x1^3 + 37/x2^3 + 19/x3^3 + 7/x4^3 + 1/x5^3 -1;
    
end