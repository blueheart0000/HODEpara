function [xo,fo] = L_SHADE(xl, options,fncobj)
% Based on
% Improving the Search Performance of SHADE Using Linear Population Size Reduction
% Ryoji Tanabe and Alex S. Fukunaga
% 2014 IEEE Congress on Evolutionary Computation (CEC) July 6-11, 2014, Beijing, China

% Version 1.01: 06/09/2024, author: Eduardo Nunes Gonçalves

% xl : vector with the initial range of the optimization variables
% opcoes: vector with options:
%         (1) initial population size
%         (2) number of objective function evaluations
%         (3) = 0 => optimization variables are free,
%             otherwise search region is options(3)*xl 

options_dflt = [100 1000 0]; % standard options 
if exist('options','var')
  options = [options, options_dflt(length(options)+1:end)];
else
  options = options_dflt;
end
  
N = options(1); 
max_NFE = options(2);
xl_mult = options(3);

Ni=N;
H=5;
Co=0.5;
Fo=0.5;
MCR=Co*ones(1,H);
MF =Fo*ones(1,H);
km =1;
Fi =zeros(N,1);
Cr =Fi;
dF =Fi;
  
n = size(xl,1); 
if exist('xi','var')
    X(:,1)=xi;
    i=2;
else
    i=1;
    xi=zeros(n,1);
end
for i=N:-1:i
  X(:,i)=xi+xl(:,1)+(xl(:,2)-xl(:,1)).*rand(n,1);
end
for i=N:-1:1
  Fx(i) = fncobj(X(:,i)');
end

NA=round(1.4*N);
Y=zeros(n,NA);
ny=0;

NFE = N; 
k=0;
[Fx, ids]=sort(Fx); 
X=X(:,ids);
fo=Fx(1);
xo = X(:,1);

u=zeros(n,1);
% 
% fprintf('%5s %5s %5s %6s %11s\n','It.','N','Na','NFE','f')
% fprintf('%5d %5d %5d %6d %11.5g\n',k,N,ny,NFE,fo) 
while NFE < max_NFE
    k=k+1;
    nb=max(2,ceil(0.11*N));
    for i = 1:N
        j=randi(H);
        Fi(i)=0;
        while Fi(i)<=0
            Fi(i)=min(1,MF(j)+0.1*tan(pi*(rand-0.5)));
        end
        % pbest index (from top nb), must differ from i
        rb = randi(nb);
        while rb == i && nb > 1
            rb = randi(nb);
        end
        % r1 from population, must differ from i
        r1 = randi(N);
        while r1 == i
            r1 = randi(N);
        end
        % r2 from population+archive, must differ from i and r1
        r2 = randi(N+ny);
        while r2 == i || r2 == r1
            r2 = randi(N+ny);
        end
        if r2<=N
            V(:,i)=X(:,i)+Fi(i)*(X(:,rb)-X(:,i)) ...
                +Fi(i)*(X(:,r1)-X(:,r2));
        else
            V(:,i)=X(:,i)+Fi(i)*(X(:,rb)-X(:,i)) ...
                +Fi(i)*(X(:,r1)-Y(:,r2-N));
        end
    end
  
    for i = 1:N
        j=randi(H);
        if MCR(j) == -1
            Cr(i) = 0;
        else
            Cr(i)=max(0,min(1,MCR(j)+0.1*randn));
        end
        di = randi(n);
        for j = 1:n
            if (j == di) || (rand <= Cr(i)) 
                u(j,1)=V(j,i);
            else
                u(j,1)=X(j,i);
            end
        end
        if xl_mult
            u = reflection(u,xl_mult*xl);
        end
        Fu = fncobj(u');
        if Fu < Fx(i)
            if ny<NA 
                ny=ny+1;
                j=ny;
            else
                j=randi(NA);
            end
            Y(:,j)=X(:,i);
            dF(i)=Fx(i)-Fu;
            X(:,i)=u;
            Fx(i)=Fu;
        else
            Cr(i)=-1;
        end
     end
    [Fx, ids]=sort(Fx);
    X=X(:,ids);
    fo=Fx(1);
    xo = X(:,1);

    NFE = NFE + N;
    Co1=0;
    Co2=0;
    Fo1=0;
    Fo2=0;
    sdF=0; 
    for i=1:N
      if Cr(i)>= 0
        sdF=sdF+dF(i);
      end
    end
    j=0;
    for i=1:N
        if Cr(i)>= 0
            j=j+1;
            di=dF(i)/sdF;
            Co1=Co1+di*Cr(i);
            Co2=Co2+di*Cr(i)^2;
            Fo1=Fo1+di*Fi(i);
            Fo2=Fo2+di*Fi(i)^2;
         end
    end
    if j>0
        if Co1 > 0
            Co=Co2/Co1;
        else
            Co=-1;
        end
        Fo=Fo2/Fo1;
    end
    MF(km)=Fo;
    MCR(km)=Co;
    km=km+1;
    if km>H
        km=1;
    end
    %fprintf('%5d %5d %5d %6d %11.5g\n',k,N,ny,NFE,fo)
    N=max(4,round((4-Ni)/max_NFE*NFE+Ni));
    NA=round(1.4*N);
    if ny > NA
       ny=NA;
    end
end

return

function xr = reflection(x,xl)
n=length(x);
xr = x;
for i=1:n
    if xr(i) < xl(i,1)
      xr(i) = xl(i,1) + (xl(i,1)-xr(i));
      if xr(i) > xl(i,2)
          xr(i)=xl(i,1);
      end
    elseif xr(i) > xl(i,2)
      xr(i) = xl(i,2) -(xr(i)- xl(i,2));
      if xr(i) < xl(i,1)
          xr(i)=xl(i,2);
      end
    end
end
return

% function f = fncobj(x)
% % Implement your objective function here
% f = 20+x(1)^2-10*(cos(2*pi*x(1)))+x(2)^2-10*(cos(2*pi*x(2))); % Rastrigin
% % Constraints:
% g = (x(1)-2)^2+(x(2)-0.5)^2-1;
% % Penalty:
% if g > 0
%     f=1e8*(1+1e8*g);
% end
% return
