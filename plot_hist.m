clear all
close all

load obj_hist_2d_CG.txt -ascii
load obj_hist_2d_MMA.txt -ascii
load obj_hist_2d_MMAPC.txt -ascii

load obj_hist_5d_CG.txt -ascii
load obj_hist_5d_MMA.txt -ascii
load obj_hist_5d_MMAPC.txt -ascii

figure(1)
set(gcf,'position',[510 180 700 400]);
plot(obj_hist_2d_CG,'-o','MarkerSize',5,'linewidth',1.5,'color',[1 0 0]);hold on;
plot(obj_hist_2d_MMA,'-o','MarkerSize',5,'linewidth',1.5,'color',[0 0.65 0]);hold on;
plot(obj_hist_2d_MMAPC,'-o','MarkerSize',5,'linewidth',1.5,'color',[0 0 1]);hold on;
set(gca,'fontsize',16,'linewidth',2.0);
xlabel('Iteration number','fontsize',16)
ylabel('{\it f}_{2D}','fontsize',16);
ylim([-1.05 -0.3])

figure(2)
set(gcf,'position',[510 180 700 400]);
plot(obj_hist_5d_CG,'-o','MarkerSize',5,'linewidth',1.5,'color',[1 0 0]);hold on;
plot(obj_hist_5d_MMA,'-o','MarkerSize',5,'linewidth',1.5,'color',[0 0.65 0]);hold on;
plot(obj_hist_5d_MMAPC,'-o','MarkerSize',5,'linewidth',1.5,'color',[0 0 1]);hold on;
set(gca,'fontsize',16,'linewidth',2.0);
xlabel('Iteration number','fontsize',16)
ylabel('{\it f}_{5D}','fontsize',16);
ylim([-1.05 -0.2])