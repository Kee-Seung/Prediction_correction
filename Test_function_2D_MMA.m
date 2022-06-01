clear all
close all

m0 = -1;
n0 = -1;
a0 = 0;
b0 = 0;

height = 0.2;
steep = -100;

m1 = height;
n1 = steep;
a1 = -0.7;
b1 = -0.1;

m2 = height;
n2 = steep;
a2 = -0.7;
b2 = 0.1;

m3 = height;
n3 = steep;
a3 = -0.5;
b3 = -0.1;

m4 = height;
n4 = steep;
a4 = -0.5;
b4 = 0.1;

objFunc = @(x) m0*exp(n0*((x(1)-a0)^2+(x(2)-b0)^2))...
    + m1*exp(n1*((x(1)-a1)^2+(x(2)-b1)^2))...
    + m2*exp(n2*((x(1)-a2)^2+(x(2)-b2)^2))...
    + m3*exp(n3*((x(1)-a3)^2+(x(2)-b3)^2))...
    + m4*exp(n4*((x(1)-a4)^2+(x(2)-b4)^2));

gradFunc = @(x) [2*m0*n0*(x(1)-a0)*exp(n0*((x(1)-a0)^2+(x(2)-b0)^2))...
    + 2*m1*n1*(x(1)-a1)*exp(n1*((x(1)-a1)^2+(x(2)-b1)^2))...
    + 2*m2*n2*(x(1)-a2)*exp(n2*((x(1)-a2)^2+(x(2)-b2)^2))...
    + 2*m3*n3*(x(1)-a3)*exp(n3*((x(1)-a3)^2+(x(2)-b3)^2))...
    + 2*m4*n4*(x(1)-a4)*exp(n4*((x(1)-a4)^2+(x(2)-b4)^2));
    2*m0*n0*(x(2)-b0)*exp(n0*((x(1)-a0)^2+(x(2)-b0)^2))...
    + 2*m1*n1*(x(2)-b1)*exp(n1*((x(1)-a1)^2+(x(2)-b1)^2))...
    + 2*m2*n2*(x(2)-b2)*exp(n2*((x(1)-a2)^2+(x(2)-b2)^2))...
    + 2*m3*n3*(x(2)-b3)*exp(n3*((x(1)-a3)^2+(x(2)-b3)^2))...
    + 2*m4*n4*(x(2)-b4)*exp(n4*((x(1)-a4)^2+(x(2)-b4)^2))];

% x_ini = [-1;0];
x_ini = [-1.01;-0.02];
x = x_ini;

warning off
UP = 3; 
LOW = -3;
rel_diff01 = 1;
rel_diff02 = 1; 
NUM_CNST = 1;
NUM_DV = length(x);
xval = x; % design variables (column vector)
xold1 = xval;
xold2 = xval;
xmin  = LOW*ones(NUM_DV,1);
xmax  = UP*ones(NUM_DV,1);
low   = xmin;
upp   = xmax;
c = 1000*ones(NUM_CNST, 1);
d = zeros(NUM_CNST, 1);
a0 = 1;
a = zeros(NUM_CNST, 1);

obj = objFunc(x);
Dobj = gradFunc(x);

f0val = obj; 
df0dx = Dobj; 
df0dx2 = 0;
   
% cnst01 = objFunc(x) + 2;
cnst01 = objFunc(x) -500;
Dcnst01 = gradFunc(x)';

% fval = [cnst01,cnst02]';
fval = cnst01';
% dfdx = [Dcnst01;Dcnst02];
dfdx = Dcnst01;
dfdx2 = 0;

loop = 0; 
MAXIT = 200;
cr01 = 10e-5;
data(1,1) = x(1);
data(1,2) = x(2);
data(1,3) = obj;
sign_list = [];
obj_hist_list = [];
while ( loop < MAXIT && (rel_diff01 > cr01 || rel_diff02 > cr01) )
    obj_hist_list(loop+1) = objFunc(x);
    f0val_old = f0val;
    loop = loop + 1;
    xval = x;

%     [sign,Z_positive,Z_negative] = pc_algorithm_2d(NUM_CNST, NUM_DV, loop, xval, xmin, xmax, xold1, xold2, low, upp, a0, a, c, d);
%     df0dx = sign*df0dx;
%     sign_list(loop) = sign;
       
      [xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, low, upp] = ...
        mmasub(NUM_CNST, NUM_DV, loop, xval, xmin, xmax, xold1, xold2, f0val, df0dx, df0dx2, fval, dfdx, dfdx2, low, upp, a0, a, c, d);
    xold2 = xold1;
    xold1 = xval;
    x = xmma; 
    
       
    %%%%%% obj cal %%%%%
    obj = objFunc(x);
    Dobj = gradFunc(x);

    f0val = obj; 
    df0dx = Dobj; 
    df0dx2 = 0;

%     cnst01 = objFunc(x) + 2;
    cnst01 = objFunc(x) -500;
    Dcnst01 = gradFunc(x)';

    fval = cnst01';
    dfdx = Dcnst01;
    dfdx2 = 0;
    
    rel_diff02 = rel_diff01;
    rel_diff01 = abs(f0val_old-f0val);

    data(loop+1,1) = x(1);
    data(loop+1,2) = x(2);
    data(loop+1,3) = obj;
    %%%%% PRINT RESULTS %%%%%
%     fprintf(['Iteration ' num2str(loop) ',  Sign : ' num2str(sign) ',  Z+ : ' num2str(Z_positive) ',  Z- : ' num2str(Z_negative) ',  f(x1,x2) = ' num2str(objFunc(x)) ',  (x1,x2) = [' num2str(x') ']\n']);
    fprintf(['Iteration ' num2str(loop) ',  f(x1,x2) = ' num2str(objFunc(x)) ',  (x1,x2) = [' num2str(x') ']\n']);
end
obj_hist_list(loop+1) = objFunc(x);

% Visualization
figure(100)
set(gcf,'renderer','painters');
n_data_x = 70;
n_data_y = 70;
x_data = linspace(-3,3,n_data_x);
y_data = linspace(-3,3,n_data_y);
obj_task2_space = zeros(n_data_x,n_data_y);
for i01 = 1:n_data_x
    for j01 = 1:n_data_y
        obj_task2_space(i01,j01) = objFunc([x_data(i01),y_data(j01)]);
    end
end
surf(x_data',y_data',obj_task2_space'); hold on;
[~,h] = contourf(x_data',y_data',obj_task2_space');hold on;
h.ContourZLevel = -4;
plot3(x_ini(1),x_ini(2),-4,'.b','markersize',15); hold on;
plot3(x(1),x(2),-4,'.r','markersize',15); hold on;
data_n = length(data);
for j01 = 1:data_n
    plot3(data(j01,1),data(j01,2),-4,'.','color',[0 0 0],'markersize',7); hold on;
    if j01 == data_n
        break
    else
        quiver3(data(j01,1),data(j01,2),-4,data(j01+1,1)-data(j01,1),data(j01+1,2)-data(j01,2),0,'r');hold on
    end
end
% xlim([-0.5 0]);
% ylim([1 1.5]);
xlim([-2 2]);
ylim([-2 2]);
zlim([-4 1]);
set(gca,'fontsize',14,'linewidth',2.0);
xlabel('{\it \chi}_1','fontsize',16);
ylabel('{\it \chi}_2','fontsize',16);
zlabel('{\it f}_{2D}','fontsize',16);
view(45,30)


figure(101)
set(gcf,'position',[510 180 400 700]);
[~,h] = contourf(x_data',y_data',obj_task2_space');hold on;
set(gca,'fontsize',16,'linewidth',4.0,'ycolor',[1 0 1],'xcolor',[1 0 1]);
h.ContourZLevel = -4;
data_n = length(data);
for j01 = 1:data_n
    plot3(data(j01,1),data(j01,2),-4,'.','color',[0 0 0],'markersize',20); hold on;
    if j01 == data_n
        break
    else
        quiver3(data(j01,1),data(j01,2),-4,data(j01+1,1)-data(j01,1),data(j01+1,2)-data(j01,2),0,'r','MaxHeadSize',0.5);hold on
    end
end
plot3(x_ini(1),x_ini(2),-4,'.b','markersize',21); hold on;
plot3(x(1),x(2),-4,'.r','markersize',21); hold on;
xlim([-1.1 0.1]);
ylim([-0.2 0.2]);
xlabel('{\it \chi}_1','fontsize',25);
ylabel('{\it \chi}_2','fontsize',25);
camroll(-90)
view([90 90])



% figure(1)
% set(gcf,'position',[510 180 700 400]);
% plot(1:loop,sign_list,'-o','MarkerSize',11,'linewidth',1.5,'color',[0 0 1]);
% set(gca,'fontsize',16,'linewidth',2.0);
% xlabel('Iteration number','fontsize',16)
% ylim([-2 3])
% xlim([1 15])
% set(gca,'ytick',[])
% xlim([0.5 15.5])
% 
figure(2)
set(gcf,'position',[510 180 700 400]);
plot(0:loop,obj_hist_list,'-o','MarkerSize',11,'linewidth',1.5,'color',[0 0 1]);
set(gca,'fontsize',16,'linewidth',2.0);
xlabel('Iteration number','fontsize',16)
ylabel('{\it f}_{2D}','fontsize',16);
% ylim([-2 3])
% xlim([1 15])

