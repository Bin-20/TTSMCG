function perf(T,logplot)

%PERF    Performace profiles
%
% PERF(T,logplot)-- produces a performace profile as described in
%   Benchmarking optimization software with performance profiles,
%   E.D. Dolan and J.J. More', 
%   Mathematical Programming, 91 (2002), 201--213.
% Each column of the matrix T defines the performance data for a solver.
% Failures on a given problem are represented by a NaN.
% The optional argument logplot is used t


%
% This function is based on the perl script of Liz Dolan.
%
% Jorge J. More', June 2004
%T=data;

if (nargin < 2) 
    logplot = 1;
end
colors  = ['g' 'r' 'c' 'b' 'm' 'k' 'y'];
%lines   = ['-.' '-' ':' ':' '-'];
linestyle = ['--' '.-''*''-' '.' '-'];
markers = ['x' '*' 's' 'd' 'v' '^' 'o'];
%np与ns分别为矩阵T的行数与列数
[np,ns] = size(T);

% Minimal performance per solver  

%解同一个问题所需时间最少的算法占用的时间，minperf是一个1*ns的矩阵。
%min(T,[],2)是一个列向量，其第i个分量为T第i行的最小值。
minperf = min(T,[],2);
% Compute ratios and divide by smallest element in each row.

r = zeros(np,ns);
for p = 1: np
  r(p,:) = T(p,:)/minperf(p); %把T中每一行的元素除以该行的最小值，得到相应的矩阵r。
end

if (logplot) 
 r = log2(r);
end


max_ratio = max(max(r)); %max_ratio为r矩阵中每一列最大值中的最大值，即矩阵r中的最大元素。

% Replace all NaN's with twice the max_ratio and
% sort，即若r中有无穷大值，则用2*max_ratio来代替。
r(find(isnan(r))) =  2*max_ratio;
%函数[r,c,v]=find(A)，其中把r与c放到一起就可以确定A中非零元的位置，v表示A的非零元。find（条件）
%表示找出满足“条件”的元素所在的位置。TF = isnan(A)返回一个与A相同维数的数组，
%若A的元素为NaN（非数值），在对应位置上返回逻辑1（真），否则返回逻辑0（假）。
r = sort(r);
%把每列从小到大排列
%B = sort(A)沿着输入参量 A的不同维的方向、从小到大重新排列 A中的元素。A 可以是字符串的、实数的
%、复数的单元数组。对于 A 中完全相同的元素，则按它们在 A 中的先后位置排列在一块；若 A 为复数的，
%则按元素幅值的从小到大排列，若有幅值相同的复数元素，则再按它们在区间[-π ,π ]的幅角从小到大排列；
%若 A 中有元素为NaN，则将它们排到最后。若 A为向量，则返回从小到大的向量，若A为二维矩阵，则按列的方向进行排列；
%若A为多维数组，sort(A)把沿着第一非单元集的元素像向量一样进行处理。


%统计好的个数及百分比
i=0;j=0;k=0;%u=0;s=0;
% l=0;u=0;
for p = 1: np
  if (r(p,1)==0.0)
    i=i+1;
  end
  if(r(p,2)==0.0)
     j=j+1;
  end    
  if(r(p,3)==0.0)
     k=k+1;
 end
% if(r(p,4)==0.0)
%   u=u+1;
% end
 %if(r(p,5)==0.0)
   %  s=s+1;
%  end
end
i
j
k
% u
%s
np
ip=i/np
jp=j/np
kp=k/np
% up=u/np
%sp=s/np
% Plot stair graphs with markers.

clf;

for s = 1: ns
     
    [xs,ys] = stairs(r(:,s),[1:np]/np);
 
% option = ['-' colors(s) markers(s)]
 
 if(s==1)
    % option=[ '-.', colors(s)]
     option=[ '-' , 'r']
 end
 if(s==2)
    % option=['-',colors(s)]
     option=['-.','b']
 end
 if(s==3)
     %   %option=[':',colors(s)]
     option=['--','m']
 end
% if(s==4)
%  %option=['',colors(s)]
%  option=['--','b']
% end
 %if(s==5)
%     % option=[':',colors(s)]
 %   option=['-','k']
%end
    %option=[colors(s),linestyle(s)]？不起作用
     % plot(xs,ys,option,'MarkerSize',3)
     
 %  xt=get(gca,'XTick');
%set(gca,'XTickLabel',sprintf('%.1f|',xt));
%set(0,'defaulttextfontsize',20);
%%%%%%下面命令中的“10”表示图形中字体大小
set(0,'defaultaxesfontsize',10);
% plot(xs,ys,option,'MarkerSize',20);
plot(xs,ys,option, 'LineWidth', 1.5)

%grid on;
 hold on;
end

% Axis properties are set so that failures are not shown,
% but with the max_ratio data points shown. This highlights
% the "flatline" effect.

%%%若是画CPU时间，应用此坐标比例
axis([ 0 1.1*max_ratio 0 1 ]);

%%%若是画迭代次数的函数时，应用此坐标比例
%axis([0  1*max_ratio  0 1 ]);
legend('TTSMCG','MPRPDP','SSDLM')
xlabel('\tau')
ylabel('\chi_S(\tau)')
%title('Fig2: Function evaluations performance profiles for the presented methods')
%title('Fig1-1: Performance profiles of with respect to CPU time in seconds')
%title('Fig3: Performance profiles with respect to the number of iterations')
%title('Fig.2-2: Performance profiles with respect to the number of iterations')
%title('Fig.2-2 Performance profiles baesd on the number of iterations')
%title('Fig.3: Performance profiles with respect to the number of function evaluations')
%title('Fig.4.1: Performance profiles with respect to the number of gradient evaluations')
%用法：axis([xmin xmax ymin ymax])
% Legends and title should be added




