


clear
load('D:\Data\Data_figure4.mat')








%% scatter ldg latency vs ncs larency
figure
idx=LDG_CSM_modulation<=0 & significant_response_any ;

scatter(LDG_latency(idx),NCS_latency(idx),12,'filled')
% title("Response latency of complex neurons")
xlabel('LDG latency (ms)')
ylabel('NCS lateny (ms)')
title([num2str(sum(idx)),' units'])
box off
ax = gca; ax.FontSize = 14; ax.FontName="Arial";
set(ax,'linewidth',1)
hl=150; h2=0:50:250;
xticks(h2);yticks(h2)
ll=min([min(NCS_latency(idx)),min(LDG_latency(idx))])-5;
xlim([-Inf,hl]); ylim([-Inf,hl])
axis square
line([ll,hl],[ll,hl],'Color','black','LineStyle','--','LineWidth',1)




%%  rectangle latency vs IGR 


IGR_1     = IGR(significant_response_any);
rect_latency_1     = rect_latency(significant_response_any);
figure
scatter(rect_latency_1,IGR_1,25,'filled')

xlabel('rectangle latency(ms)'); ylabel('IGR')
h1=lsline; h1.LineWidth = 2;
ax = gca; ax.FontSize = 14; ax.FontName="Arial";

%% scatter surround modulation vs IGR
figure


scatter(surround_modulation(significant_response_any),IGR_1,15,'filled','MarkerFaceColor',[ 0 0.4470 0.7410])
ylabel('IGR'); xlabel('surround modulation')

ax = gca; ax.FontSize = 18; ax.FontName="Arial"; set(ax,'linewidth',1.5)
axis square; h1 = refline; h1.LineWidth = 2;



%% histogram for surround_suppression and separated igr
idx=IGR<0 & significant_response_any;
idx2=IGR>=0 & significant_response_any;
ldg_pref_med=nanmedian(surround_modulation(idx));
ncs_pref_med=nanmedian(surround_modulation(idx2));



figure;hold on
histogram(surround_modulation(idx),25,'Normalization','probability','FaceColor','red');
histogram(surround_modulation(idx2),25,'Normalization','probability','FaceColor','blue');
xli1=line([ldg_pref_med ldg_pref_med],[0 .14],'Color','r','LineStyle','--','LineWidth',1.5);
text1=text(ldg_pref_med,-.9,'\leftarrow median');text1.FontSize = 16;
xli2=line([ncs_pref_med ncs_pref_med],[0 .14],'Color','b','LineStyle','--','LineWidth',1.5);
text2=text(ncs_pref_med,.9,'\leftarrow median'); text2.FontSize = 16;
set(get(get(xli1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(xli2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
xlabel('surround modulation'); ylabel('probability');box off;
legend('LDG-preferring','NCS-preferring','Location','northwest')
ax = gca; ax.FontSize = 18; ax.FontName="Arial";set(ax,'linewidth',1.5)



%% IGRF1 vs CSM

x  = LDG_CSM_modulation(significant_response_both);
y  = IGRf1(significant_response_both);


figure;
scatter(x,y,15,'filled')
ylabel('IGR_{F1}');xlabel('CSM index');box off
h1 = refline; h1.LineWidth = 2;
ax = gca; ax.FontSize = 14; ax.FontName="Arial";axis square


%% histogram for TPlatency

figure%('Visible','off');
histogram(TPlatency(significant_response_any),30);
xlabel('TP latency (ms)')
ylabel('# units')
% grid on
ax = gca;
ax.FontSize = 18; 
ax.FontName="Arial";
xlim([0.126,1.182])


%% Icell/Ecell
X = categorical({'I cells','E cells'});
X = reordercats(X,{'I cells','E cells'});
%% I/E percentage

Icell_ncs_evoked=sum(Icell_ind & significant_response_NCS)/sum(Icell_ind)*100;
Icell_ldg_evoked=sum(Icell_ind & significant_response_LDG)/sum(Icell_ind)*100;
Ecell_ncs_evoked=sum(Ecell_ind & significant_response_NCS)/sum(Ecell_ind)*100;
Ecell_ldg_evoked=sum(Ecell_ind & significant_response_LDG)/sum(Ecell_ind)*100;
barr=[Icell_ncs_evoked,Icell_ldg_evoked;Ecell_ncs_evoked,Ecell_ldg_evoked];
figure
bar(X,barr)
ylabel('% of evoked cells')
legend('NCS evoked','LDG evoked','location','northwestoutside')
ax = gca; ax.FontSize = 14;  ax.FontName="Arial";set(ax,'linewidth',1)


%% IE vs IGRf1
IGRf1_Icell=IGRf1(Icell_ind);
IGRf1_Ecell=IGRf1(Ecell_ind);

pd_Icell = fitdist(IGRf1_Icell,'Normal');
pd_Ecell = fitdist(IGRf1_Ecell,'Normal');
ci_Icell = pd_Icell.paramci;
ci_Ecell = pd_Ecell.paramci;
ci_mean  = [ci_Icell(:,1)-pd_Icell.mean,ci_Ecell(:,1)-pd_Ecell.mean];
Y       = [pd_Icell.mean,pd_Ecell.mean];
% [~,p2]=ttest2(IGRf1_Ecell,IGRf1_Icell);
figure;hold on
bar(X,Y,0.5)
errorbar(X,Y,ci_mean(1,:),ci_mean(2,:),'LineWidth', 2,'CapSize',20,'LineStyle', 'none')
ylabel('IGR_{F1}'); box on
ax = gca;ax.FontSize = 18; ax.FontName="Arial";

%% IE VS igr
IGR_Icell=IGR(Icell_ind);
IGR_Ecell=IGR(Ecell_ind);

pd_Icell = fitdist(IGR_Icell,'Normal');
pd_Ecell = fitdist(IGR_Ecell,'Normal');
ci_Icell = pd_Icell.paramci;
ci_Ecell = pd_Ecell.paramci;
ci_mean  = [ci_Icell(:,1)-pd_Icell.mean,ci_Ecell(:,1)-pd_Ecell.mean];
Y       = [pd_Icell.mean,pd_Ecell.mean];
% [~,p2]=ttest2(IGR_Ecell,IGR_Icell)
figure;hold on
bar(X,Y,0.5)
errorbar(X,Y,ci_mean(1,:),ci_mean(2,:),'LineWidth', 2,'CapSize',20,'LineStyle', 'none')
ylabel('IGR')


%% full screen experiment results

% IE VS igr non overlapping
fullfeild_IGR_Ecell=fullfeild_IGR(fullfeild_Ecell_ind & fullfeild_non_overlapping_V1 & fullfeild_significant_response_any);
fullfeild_IGR_Icell=fullfeild_IGR(~fullfeild_Ecell_ind & fullfeild_non_overlapping_V1& fullfeild_significant_response_any);

pd_Icell = fitdist(fullfeild_IGR_Icell,'Normal');
pd_Ecell = fitdist(fullfeild_IGR_Ecell,'Normal');
ci_Icell = pd_Icell.paramci;
ci_Ecell = pd_Ecell.paramci;
ci_mean  = [ci_Icell(:,1)-pd_Icell.mean,ci_Ecell(:,1)-pd_Ecell.mean];
Y       = [pd_Icell.mean,pd_Ecell.mean];
figure;hold on
bar(X,Y,0.5)
errorbar(X,Y,ci_mean(1,:),ci_mean(2,:),'LineWidth', 2,'CapSize',20,'LineStyle', 'none')
ylabel('IGR')
ax = gca; ax.FontSize = 18; ax.FontName="Arial";


% IE vs IGRf1 non overlapping
fullfeild_IGRf1_Ecell=fullfeild_IGRf1(fullfeild_Ecell_ind & fullfeild_non_overlapping_V1& fullfeild_significant_response_any);
fullfeild_IGRf1_Icell=fullfeild_IGRf1(~fullfeild_Ecell_ind & fullfeild_non_overlapping_V1& fullfeild_significant_response_any);
pd_Icell = fitdist(fullfeild_IGRf1_Icell,'Normal');
pd_Ecell = fitdist(fullfeild_IGRf1_Ecell,'Normal');
ci_Icell = pd_Icell.paramci;
ci_Ecell = pd_Ecell.paramci;
ci_mean  = [ci_Icell(:,1)-pd_Icell.mean,ci_Ecell(:,1)-pd_Ecell.mean];
Y       = [pd_Icell.mean,pd_Ecell.mean];
figure;hold on
bar(X,Y,0.5)
errorbar(X,Y,ci_mean(1,:),ci_mean(2,:),'LineWidth', 2,'CapSize',20,'LineStyle', 'none')
ylabel('IGR_{F1}'); box on
ax = gca;ax.FontSize = 18; ax.FontName="Arial";



%% I-E delay VS RF dist 
idx1 =fullfeild_significant_response_any  & fullfeild_NCS_CSM_modulation<0 &fullfeild_Icell_ind;
idx2 =fullfeild_significant_response_any  & fullfeild_NCS_CSM_modulation<0 & ~fullfeild_Icell_ind;
figure
x = fullfeild_RF_distance(idx1); y  = fullfeild_ncs_delay(idx1);
x= 2*atand(x./24); % convert to deg
scatter(x,y,25,"diamond",'filled'); hold on;

refline
xlabel('RF/inducer distance (degrees)');ylabel('NCS delay index')

ax = gca; ax.FontSize = 12; ax.FontName="Arial";
legend(['#',num2str(sum(~isnan(fullfeild_ncs_delay(idx1)))),' I cell'],"Box","off",'Location','northeastoutside','FontSize',13);axis square
ylim([-1,1])

figure

x = fullfeild_RF_distance(idx2); y  = fullfeild_ncs_delay(idx2);
x= 2*atand(x./24); % convert to deg

scatter(x,y,25,"hexagram",'filled')
refline;
xlabel('RF/inducer distance (degree)');ylabel('NCS delay index')

ax = gca; ax.FontSize = 12; ax.FontName="Arial";
legend(['#',num2str(sum(~isnan(fullfeild_ncs_delay(idx2)))),' E cell' ...
    ],"Box","off",'Location','northeastoutside','FontSize',13);axis square
ylim([-1,1])  % xlim([0,10])
