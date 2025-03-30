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
%np��ns�ֱ�Ϊ����T������������
[np,ns] = size(T);

% Minimal performance per solver  

%��ͬһ����������ʱ�����ٵ��㷨ռ�õ�ʱ�䣬minperf��һ��1*ns�ľ���
%min(T,[],2)��һ�������������i������ΪT��i�е���Сֵ��
minperf = min(T,[],2);
% Compute ratios and divide by smallest element in each row.

r = zeros(np,ns);
for p = 1: np
  r(p,:) = T(p,:)/minperf(p); %��T��ÿһ�е�Ԫ�س��Ը��е���Сֵ���õ���Ӧ�ľ���r��
end

if (logplot) 
 r = log2(r);
end


max_ratio = max(max(r)); %max_ratioΪr������ÿһ�����ֵ�е����ֵ��������r�е����Ԫ�ء�

% Replace all NaN's with twice the max_ratio and
% sort������r���������ֵ������2*max_ratio�����档
r(find(isnan(r))) =  2*max_ratio;
%����[r,c,v]=find(A)�����а�r��c�ŵ�һ��Ϳ���ȷ��A�з���Ԫ��λ�ã�v��ʾA�ķ���Ԫ��find��������
%��ʾ�ҳ����㡰��������Ԫ�����ڵ�λ�á�TF = isnan(A)����һ����A��ͬά�������飬
%��A��Ԫ��ΪNaN������ֵ�����ڶ�Ӧλ���Ϸ����߼�1���棩�����򷵻��߼�0���٣���
r = sort(r);
%��ÿ�д�С��������
%B = sort(A)����������� A�Ĳ�ͬά�ķ��򡢴�С������������ A�е�Ԫ�ء�A �������ַ����ġ�ʵ����
%�������ĵ�Ԫ���顣���� A ����ȫ��ͬ��Ԫ�أ��������� A �е��Ⱥ�λ��������һ�飻�� A Ϊ�����ģ�
%��Ԫ�ط�ֵ�Ĵ�С�������У����з�ֵ��ͬ�ĸ���Ԫ�أ����ٰ�����������[-�� ,�� ]�ķ��Ǵ�С�������У�
%�� A ����Ԫ��ΪNaN���������ŵ������ AΪ�������򷵻ش�С�������������AΪ��ά�������еķ���������У�
%��AΪ��ά���飬sort(A)�����ŵ�һ�ǵ�Ԫ����Ԫ��������һ�����д���


%ͳ�ƺõĸ������ٷֱ�
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
    %option=[colors(s),linestyle(s)]����������
     % plot(xs,ys,option,'MarkerSize',3)
     
 %  xt=get(gca,'XTick');
%set(gca,'XTickLabel',sprintf('%.1f|',xt));
%set(0,'defaulttextfontsize',20);
%%%%%%���������еġ�10����ʾͼ���������С
set(0,'defaultaxesfontsize',10);
% plot(xs,ys,option,'MarkerSize',20);
plot(xs,ys,option, 'LineWidth', 1.5)

%grid on;
 hold on;
end

% Axis properties are set so that failures are not shown,
% but with the max_ratio data points shown. This highlights
% the "flatline" effect.

%%%���ǻ�CPUʱ�䣬Ӧ�ô��������
axis([ 0 1.1*max_ratio 0 1 ]);

%%%���ǻ����������ĺ���ʱ��Ӧ�ô��������
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
%�÷���axis([xmin xmax ymin ymax])
% Legends and title should be added




