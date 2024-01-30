

%% generating figure 2
addpath(genpath('D:\matwork\circstat'))

load('D:\Data\NCS_2024_figure2.mat')

Cond_names  = { 'NCS';'LDG';'DBC'};


%% scatter ldg vs ncs
figure
scatter(LDG_FR(significant_response_any),NCS_FR(significant_response_any),12,'filled')
xlabel('LDG');ylabel('NCS')
ax = gca; ax.FontSize = 18; ax.FontName="Arial"; set(ax,'linewidth',1.5)
xlim([-Inf,60]); ylim([-Inf,60])
x = [min([min(LDG_FR(significant_response_any)),min(NCS_FR(significant_response_any))]),60];
line(x,x,'Color','black','LineStyle','--','LineWidth',1)
xline(0,'LineWidth',1); yline(0,'LineWidth',1);
axis square



%% scatter ldg vs dbc
figure
scatter(LDG_FR(significant_response_any),DBC_FR(significant_response_any),12,'filled')
xlabel('LDG');ylabel('DBC')
ax = gca; ax.FontSize = 18; ax.FontName="Arial"; set(ax,'linewidth',1.5)
xlim([-Inf,60]); ylim([-Inf,60])
x = [min([min(LDG_FR(significant_response_any)),min(DBC_FR(significant_response_any))]),60];
line(x,x,'Color','black','LineStyle','--','LineWidth',1)
xline(0,'LineWidth',1); yline(0,'LineWidth',1);
axis square



%% barplot firing rate

%%% first we do perform hierarchical bootstrapping
parm = 'mean';
nboot = 10000;

data = [mouse_day_idx(significant_response_any,1),NCS_FR(significant_response_any)];
ncs_bootstat = f_bootstrapped_samples(data,nboot,parm);
data = [mouse_day_idx(significant_response_any,1),LDG_FR(significant_response_any)];
ldg_bootstat = f_bootstrapped_samples(data,nboot,parm);
data = [mouse_day_idx(significant_response_any,1),DBC_FR(significant_response_any)];
dbc_bootstat = f_bootstrapped_samples(data,nboot,parm);

b=[ncs_bootstat,ldg_bootstat,dbc_bootstat];
sem_boot = std(b); 
mean_boot= mean(b);

[~,~,stats] = kruskalwallis(b);

X = categorical({'NCS','LDG', 'DBC'});
X = reordercats(X,{'NCS','LDG', 'DBC'});

figure('units','normalized','outerposition',[0.3 0.3 0.3 0.7])
subplot(2,1,1);hold on
bar(X,mean_boot,0.5)
errorbar(X,mean_boot,sem_boot,'black','LineWidth', 2,'CapSize',20,'LineStyle', 'none')
ylabel('Spike/sec')
ax = gca; ax.FontSize = 14; ax.FontName="Arial"; set(ax,'linewidth',1.5)
box off; axis square

c = multcompare(stats,'Display','off');
subplot(2,1,2);hold on
scatter(c(c(:,6)<0.05,1),c(c(:,6)<0.05,2),505,'hb','filled')
scatter(c(c(:,6)<0.05,2),c(c(:,6)<0.05,1),505,'hb','filled')
ylim([0.5,3.5]);yticks([1,2,3]);yticklabels(Cond_names)
xlim([0.5,3.5]);xticks([1,2,3]);xticklabels(Cond_names)
title({'post-hoc multiple comparison test'})
grid on;axis square;ax = gca;ax.FontSize = 14;ax.FontName = "Arial";



%% histogram for angle diff
figure;
pref_chang = rad2deg(wrapToPi(deg2rad(pref_ang_LDG(significant_response_both)-pref_ang_NCS(significant_response_both))));
histogram(pref_chang,-186:15:185)
% title(['changes in preferred angle for ',area_names{area}])
% xlabel("LDG_{preferred angle}-NIG_{preferred angle}")
xlabel("changes in preferred angle")
ylabel("# neurons")
ax=gca;
ax.FontSize = 18;
ax.FontName="Arial";
set(ax,'linewidth',1.5)
box off
axis square
% ylim([0,450])
xticks(-180:90:180)




%% polar histogram for phase diference

b=deg2rad(Phase_diff(F1_dominants));
% circ_R= circ_r(b);
% [pval, z] = circ_rtest(b)
[mu, ul, ll] = circ_mean(b);
% [s, s0] = circ_std(b);

figure
polhist=polarhistogram(b,71,'EdgeColor','blue','LineWidth',0.2);
hold on
polarplot([0,mu],[0,max(polhist.Values)],'r','LineWidth',3)
title([num2str(sum(F1_dominants)),' units'])

ax = gca;
ax.FontSize = 14; 
ax.FontName="Arial";
set(ax,'linewidth',1)



