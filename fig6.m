clear
load('D:\Data\Data_figure6.mat')


diameter_change=[pup_dilation_index(NCS_trials_p), ones(sum(NCS_trials_p),1);...
                 pup_dilation_index(LDG_trials_p), ones(sum(LDG_trials_p),1)+1;...
                 pup_dilation_index(DBC_trials_p), ones(sum(DBC_trials_p),1)+2];
     
mean_diameter_change =[nanmean(pup_dilation_index(NCS_trials_p)),...
                       nanmean(pup_dilation_index(LDG_trials_p)),...
                       nanmean(pup_dilation_index(DBC_trials_p))];
         
std_diameter_change = [std(pup_dilation_index(NCS_trials_p),'omitnan')/sqrt(sum(NCS_trials_p)),...
                       std(pup_dilation_index(LDG_trials_p),'omitnan')/sqrt(sum(LDG_trials_p)),...
                       std(pup_dilation_index(DBC_trials_p),'omitnan')/sqrt(sum(DBC_trials_p))];


[~,~,stats] = kruskalwallis(diameter_change(:,1),diameter_change(:,2),'off');

[c,m] = multcompare(stats,'Display','off');

figure('units','normalized','outerposition',[0 0 .2 .7])
subplot(2,1,1)
bar(mean_diameter_change)
hold on
errorbar(mean_diameter_change,std_diameter_change,'b','CapSize',0,'LineWidth',1,'CapSize',20,'LineStyle', 'none');

xticklabels({'NCS';'LDG';'DBC'})
ylabel('Pupil Dilation Index')
ax = gca; ax.FontSize = 12;

subplot(2,1,2)
scatter(c(c(:,6)<0.05,1),c(c(:,6)<0.05,2),200,'hb','filled')
hold on;box on
scatter(c(c(:,6)<0.05,2),c(c(:,6)<0.05,1),200,'hb','filled')
ylim([0.5,3.5]);yticks([1,2,3]);yticklabels({'NCS';'LDG';'DBC'})
xlim([0.5,3.5]);xticks([1,2,3]);xticklabels({'NCS';'LDG';'DBC'})
% title('ranksum test')
title('post-hoc multi comparison' )

grid on;axis square;  ax = gca; ax.FontSize = 12;ax.FontName = "Arial";



%% normalize pupil size
Pupil_time_ax(1,:) = (Pupil_time_ax(2,:)) + 5*diff(Pupil_time_ax(1,1:2));
pre_med   = (nanmedian(pupil_diameter_f(:,1:5),2));
pupil_diameter_f2 = (pupil_diameter_f-pre_med)./pre_med.*100;%./max_pupil;

% pre_pup_size  = nanmean(diameter_data_f(:,5));
% pupil_diameter_f2 = pupil_diameter_f2 - pre_pup_size;




trial_discard_idx = sum(isoutlier(pupil_diameter_f,2),2)>0 | sum(isoutlier(pupil_diameter_f),2)>0 ;%& diameter_data_f(:,6);
ncs_idx = NCS_trials_p & ~trial_discard_idx;
ldg_idx = LDG_trials_p & ~trial_discard_idx;
dbc_idx = DBC_trials_p & ~trial_discard_idx;
tr_idx  = [ncs_idx,ldg_idx,dbc_idx];

% & dm_change > 0 & dm_change < .2;
% ncs_dia = (pupil_diameter_f(ncs_idx,:));
% ldg_dia = (pupil_diameter_f(ldg_idx,:));
% dbc_dia = (pupil_diameter_f(dbc_idx,:));

% pre_med   = nanmean(nanmedian(ncs_dia(:,1:7),2));
% ncs_dia = (ncs_dia- pre_med)./pre_med.*100; 
% 
% pre_med   = nanmean(nanmedian(ldg_dia(:,1:7),2));
% ldg_dia = (ldg_dia- pre_med)./pre_med.*100;
% 
% pre_med   = nanmean(nanmedian(dbc_dia(:,1:7),2));
% dbc_dia = (dbc_dia- pre_med)./pre_med.*100;

mean_ncs_dia = nanmean(pupil_diameter_f2(ncs_idx,:));%-nanmedian(nanmedian(pupil_diameter_f2(ncs_idx,1:15)));
mean_ldg_dia = nanmean(pupil_diameter_f2(ldg_idx,:));%-nanmedian(nanmedian(pupil_diameter_f2(ldg_idx,1:15)));
mean_dbc_dia = nanmean(pupil_diameter_f2(dbc_idx,:));%-nanmedian(nanmedian(pupil_diameter_f2(dbc_idx,1:15)));


% mean_ncs_dia = nanmean(ncs_dia);
% mean_ldg_dia = nanmean(ldg_dia);
% mean_dbc_dia = nanmean(dbc_dia);
mean_diameter= [mean_ncs_dia;mean_ldg_dia;mean_dbc_dia];

std_ncs_dia = std(pupil_diameter_f2(ncs_idx,:),'omitnan')/sqrt(sum(tr_idx(:,1)));
std_ldg_dia = std(pupil_diameter_f2(ldg_idx,:),'omitnan')/sqrt(sum(tr_idx(:,2)));
std_dbc_dia = std(pupil_diameter_f2(dbc_idx,:),'omitnan')/sqrt(sum(tr_idx(:,3)));
std_diameter= [std_ncs_dia;std_ldg_dia;std_dbc_dia];

%% compute confidence interval
ci_diameter = nan(3,length(Pupil_time_ax(1,:)));
ci_diameterl = nan(2,length(Pupil_time_ax(1,:)));
ci_mean=zeros(3);
ci_meanl=zeros(2);


for i=1:length(Pupil_time_ax(1,:))
    
    pd = fitdist(pupil_diameter_f2(ncs_idx,i),'Normal');
    ci = pd.paramci-pd.mean;
    ci_mean(1,i) =ci(2,1);
    pd = fitdist(pupil_diameter_f2(ldg_idx,i),'Normal');
    ci = pd.paramci-pd.mean;
    ci_mean(2,i) =ci(2,1);
    pd = fitdist(pupil_diameter_f2(dbc_idx,i),'Normal');
    ci = pd.paramci-pd.mean;
    ci_mean(3,i) =ci(2,1);
    


end

alph=.5;
 
figure
h1 = errorbar(Pupil_time_ax(1,1:end),mean_diameter(1,1:end),ci_mean(1,1:end),'b','CapSize',3,'LineWidth',1.5);
set(h1.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [h1.Line.ColorData(1:3); 255*alph])
hold on
h2 = errorbar(Pupil_time_ax(1,1:end),mean_diameter(2,1:end),ci_mean(2,1:end),'r','CapSize',3,'LineWidth',1.5);
set(h2.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [h2.Line.ColorData(1:3); 255*alph])
h3 = errorbar(Pupil_time_ax(1,1:end),mean_diameter(3,1:end),ci_mean(3,1:end),'m','CapSize',3,'LineWidth',1.5);
set(h3.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [h3.Line.ColorData(1:3); 255*alph])
legend('NCS','LDG','DBC','Location','southeast')
legend('boxoff')

xl = xline(0,'LineWidth',1,'Color','g');
set(get(get(xl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
xl = xline(1,'LineWidth',1,'Color','g');
set(get(get(xl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
xlabel('time(s)')
ylabel('pupil diameter change (%)')
ax = gca; ax.FontSize = 12;ax.FontName = "Arial";

