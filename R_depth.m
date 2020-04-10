% 刻蚀深度反射曲线
% 使用说明：
% 1， 在“用户编写部分”设置材料层状结构；
% 2，运行此程序，得到反射率与厚度的关系
% 作者：谭少阳 时间：2015-12-21  E-mail：tanshy10@semi.ac.cn
clc;
clear;
load meterial;
%%
lambda=678; % nm 监测激光波长

N_air=1;


%% 材料层状整体结构

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%用户编写部分%%%%%%%%%%%%%%
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% 折射率参数
N_gaas=interp1(gaas(:,1),gaas(:,2)-1i*gaas(:,3),lambda);
N_algaas1=3.432-0.00168i;   % 0.47
N_algaas2=3.647-0.0706i;    %0.26
N_algaas3=3.731-0.1049i;    % 0.1
N_ingaas=interp1(ingaas(:,1),ingaas(:,2)-1i*ingaas(:,3),lambda);
% 材料层状结构
% 注意：h0是包含了空气与沉底的整体结构，其中空气和沉底的厚度均为0；
h0=[0,      300,   500,       500,         50,50               800,        50,     0]; % 模拟h_l,h_h;    
n=[N_air,   N_gaas,N_algaas2,  N_algaas3,  N_algaas2,N_gaas,   N_algaas3,  N_gaas, N_ingaas];
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %%%%%%%%%%%%用户编写部分%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_angle=degtorad(0); % 测量系统入射角
Depth=[];
r_E=[];

for i_l=2:length(h0)-1
    h=h0;
    h(2:i_l-1)=0; 
    
% 1、计算反射曲线
x_var=0:(h0(i_l)-1);    % 分辨率为1 nm
fun_var=0*x_var;

for i_var=1:length(x_var)
h(i_l) = h0(i_l)-x_var(i_var);
theta=asin(sin(i_angle)*n(1)./n);        %光在各层内的角度
delta=2*pi/lambda*n.*h.*cos(theta);
eta=n.*cos(theta);
r=zeros(length(n)-1);
for i=1:length(n)-1
    r(i)=(eta(i)-eta(i+1))/(eta(i)+eta(i+1));
end
mu=1;
for i=1:length(r);
    m=[exp(1i*delta(i)),r(i)*exp(1i*delta(i));r(i)*exp(-1i*delta(i)),exp(-1i*delta(i))];
    mu=mu*m;
end
fun_var(i_var)=mu(2,1)/mu(1,1);
end

Depth=[Depth,sum(h0(1:i_l-1))+x_var];
r_E=[r_E,fun_var];
end
%
R=abs(r_E).^2;
phase=angle(r_E)*180/pi;
dR=10*diff(R,1);
dR=[dR(1),dR];

% 2、绘制反射曲线
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultAxesLineWidth',2);
figure1 = figure;
axes1 = axes('Parent',figure1,'XMinorTick','on');
box(axes1,'on');
hold(axes1,'all');
plot(Depth,R,'r-','linewidth',2);

ylabel('reflection amplitude','fontsize',20);
xlabel('depth (nm)','fontsize',20);
% axis([0 5000 0.25 0.45]);

% 绘制界面标示
for i=1:length(h0)
    temp=sum(h0(1:i));
    plot([temp,temp],[0.3,0.36],'k-');
end

% subplot(1,2,2)
% plot(lambda,phase);
% title('reflection phase');
% xlabel('wavelength (\mu{m})');
