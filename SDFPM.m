function SDFPM( )
% this master file uses cg method of Hager & Zhang to solve nonlinear
% monotone equations with convex constraints
%
% This version Dec. 2, 2013, by Jinkui Liu

%%%%%%%%%%%%%%%%%问题选取及维数设置%%%%%%%%%%%
tic;
problem =  2;% =1,2,3,4,5
n =100000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%获取初始点%%%%%%%%%%%%%%%%
%%%%%obtain the given initial point whose each component is different

b=ones(n,1);
for i=1:n
b(i)=0.1;
end
xk=b;

%%%%%obtain the initial point by random selection 
% xk =-1+2*rand(n,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%参数初始化%%%%%%%%%%%%%%%
rho = 0.55;
sigma = 10^(-6);
beta=1;
iter = 0;
nf = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%算法开始执行%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%调用测试函数%%%%%%%%%%%%%%
if problem==1
    F=@F_P1;
    Pc=@Pc_P1;
elseif problem==2
    F=@F_P2;
    Pc=@Pc_P2;
elseif problem==3
    F=@F_P3;
    Pc=@Pc_P3; 
elseif problem==4
    F=@F_P4;
    Pc=@Pc_P4; 
elseif problem==5
    F=@F_P5;
    Pc=@Pc_P5; 
elseif problem==6
    F=@F_P6;
    Pc=@Pc_P6;
elseif problem==7
    F=@F_P7;
    Pc=@Pc_P7;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fk = F(xk);
nf = nf+1;
dk = -fk;
%%%%%%%%%%%%%%%%%开始迭代%%%%%%%%%%%%%%%%%%
while norm(fk)>10^(-5)
    tk=beta;
    zk = xk + tk*dk;
    fzk = F(zk);
    nf=nf+1;
    while (-fzk'*dk < sigma*tk*norm(fzk)*norm(dk)^2 & tk>10^(-10))
        tk = rho*tk;
        zk = xk + tk*dk;
        fzk = F(zk);
        nf = nf+1;      
    end
%     sk=tk*dk;
    alphak = fzk'*(xk-zk)/norm(fzk)^2;
    xk1 = xk - 1.2*alphak*fzk;
    xknew = Pc(xk1);
    fknew = F(xknew);
    nf = nf+1;
    yk = fknew-fk;
    sk=xknew-xk;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ck=norm(sk)/norm(fk)+max(0,-(fk'*yk)/(fk'*fk));
vk=yk+ck*fk;
mk=max(sk'*yk,max(dk'*dk,fk'*vk));
ak=(fknew'*yk)/mk-(2*fknew'*sk)/mk+((fknew'*sk)*(sk'*fk))/(mk*(sk'*sk));
bk=(fknew'*sk)/mk-((fknew'*sk)*(sk'*fk)*(sk'*yk))/((sk'*sk)*(mk)^2);

if (fknew'*sk<=0) 
    if (fknew'*fk>=0)
    dknew=-fknew+ak*sk+bk*fk;
    else
    dknew=-fknew+ak*sk;
    end
 else
    dknew=-fknew;
end
    %%%%Preparations for the next iteration
    fk = fknew;   
    xk = xknew;
    dk = dknew;
    iter = iter +1; 
    
end
toc;
%output file  
diary datamycg_descent_eq.m
      sprintf('%5i&%5i& %5i& %.4f&  %.5e\n',n,iter,nf,toc,norm(fk))     
      diary off
%%%%%%%%%%%%%%%%%%%%%%算法终止%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
%%---------------Define F(*) and Pc(*)--------------------


%------Problem 1------
function y=F_P1(x)
n=length(x);
y=ones(n,1);
for i=1:n
    y(i)=(exp(x(i)))^2+3*sin(x(i))*cos(x(i))-1;
 end   
%Projection operator
function y=Pc_P1(x)
y=max(x,-5);
%------Problem 2------
function y=F_P2(x)
n=length(x);
y=ones(n,1);
y(1)=4*x(1)-x(2)+exp(x(1))-5;
y(n)=-x(n-1)+4*x(n)+exp(x(n))-5;
for ii=2:(n-1)
    y(ii)=-x(ii-1)+4*x(ii)-x(ii+1)+exp(x(ii))-5;
 end      
%Projection operator
function y=Pc_P2(x)
%n=length(x);
%if sum(x)<=n && min(x)>=-1
 %   y=x;
%else
%    y=quadprog(eye(n),-x,[ones(1,n);zeros(n-1,n)],[n;zeros(n-1,1)],[],[],-ones(n,1),[]);
%end
y=max(x,-1);
%------Problem 3------
function y=F_P3(x)
n=length(x);
y=ones(n,1);
h=1/(n+1);
         y(1)=x(1)-exp(cos((x(1)+x(2))/1));
         y(n)=x(n)-exp(cos(h*(x(n-1)+x(n))/n));
         for i=2:(n-1)
             y(i)=x(i)-exp(cos(h*(x(i-1)+x(i)+x(i+1))/i));
         end           
%Projection operator
function y=Pc_P3(x)
y=max(x,0);


%------Problem 4------
function y=F_P4(x)

n=length(x);
y=ones(n,1);
         h=1/(n+1);
         y(1)=2*x(1)+0.5*h^2*(x(1)+h)^3-x(2);
         y(n)=2*x(n)+0.5*h^2*(x(n)+h*n)^3-x(n-1);
         for i=2:(n-1)
            y(i)=2*x(i)+0.5*h^2*(x(i)+h*i)^3-x(i-1)+x(i+1);
         end   

%Projection operator
function y=Pc_P4(x)
    
y=max(x,-1);


%------Problem 5------
function y=F_P5(x)
%Monotone nonlinear equations
n=length(x);
y=ones(n,1);
y(n)=x(n)-0.1*x(1)^2;
for i=1:(n-1)
    y(i)=x(i)-0.1*x(i+1)^2;
end
%Projection operator
function y=Pc_P5(x)   
y=max(x,-5);

%------Problem 6------
 function y=F_P6(x)
%Monotone nonlinear equations
n=length(x);
y=ones(n,1);
y(1)=4*x(1)*(x(1)^2+x(2)^2)-5;
for i=2:n-1
    y(i)=4*x(i)*(x(i-1)^2+x(i)^2)+4*x(i)*(x(i)^2+x(i+1)^2)-5;
end
y(n)=4*x(n)*(x(n-1)^2+x(n)^2)-5;
%Projection operator
function y=Pc_P6(x)
    
y=max(x,0);

 
%------Problem 7------
function y=F_P7(x)
n=length(x);
y=ones(n,1);
for i=1:n
y(i)=x(i)-3*x(i)*(sin(x(i))/3-0.66)+2;
end
%Projection operator
function y=Pc_P7(x)
 
y=max(x,-5);     