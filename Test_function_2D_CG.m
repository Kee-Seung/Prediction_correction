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

obj = @(x) m0*exp(n0*((x(1)-a0)^2+(x(2)-b0)^2))...
    + m1*exp(n1*((x(1)-a1)^2+(x(2)-b1)^2))...
    + m2*exp(n2*((x(1)-a2)^2+(x(2)-b2)^2))...
    + m3*exp(n3*((x(1)-a3)^2+(x(2)-b3)^2))...
    + m4*exp(n4*((x(1)-a4)^2+(x(2)-b4)^2));

Dobj = @(x) [2*m0*n0*(x(1)-a0)*exp(n0*((x(1)-a0)^2+(x(2)-b0)^2))...
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

f0val = obj(x);
f0val_old = f0val + 1;

df0dx  = Dobj(x);
dir = -df0dx;

loop = 0;
MAXIT = 500;
cr01 = 1e-5;

data(1,1) = x(1);
data(1,2) = x(2);
data(1,3) = obj(x);
alpha_old = 0.0025;
dir_old = dir;
Dobj_old = df0dx;
obj_hist_list = [];
while loop < MAXIT && abs((f0val_old-f0val)/f0val)>cr01 && norm(df0dx)>cr01
    obj_hist_list(loop+1) = obj(x);
    
    % iteration counter
    loop = loop + 1;
    
    % Step length
    alpha = Backtracking_Line_Search(obj,f0val,x,df0dx,Dobj_old,dir,dir_old,alpha_old);
            
    % update x
    x = x + alpha*dir;
    alpha_old = alpha;

    % update obj func values
    f0val_old = f0val;
    f0val = obj(x);
    
    % update dx
    Dobj_old = df0dx;
    df0dx    = Dobj(x);
    
    % Beta value in Fletcher-Reeves method
    beta = (df0dx'*df0dx)/(Dobj_old'*Dobj_old);
    
    % update search direction
    dir_old = dir;
    dir = -df0dx + beta*dir;
    
    data(loop+1,1) = x(1);
    data(loop+1,2) = x(2);
    data(loop+1,3) = obj(x);
    
    % Print output on screen
    history = (['Iteration: ' sprintf('%d',loop) '   Obj: ' sprintf('%1.4f',f0val) '   x: ' sprintf('% 1.4f',x')]);
    disp(history)
end
obj_hist_list(loop+1) = obj(x);

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
        obj_task2_space(i01,j01) = obj([x_data(i01),y_data(j01)]);
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


figure(2)
set(gcf,'position',[510 180 700 400]);
plot(0:loop,obj_hist_list,'-o','MarkerSize',6,'linewidth',1.5,'color',[0 0 1]);
set(gca,'fontsize',16,'linewidth',2.0);
xlabel('Iteration number','fontsize',16)
ylabel('{\it f}_{2D}','fontsize',16);
xlim([0 100])