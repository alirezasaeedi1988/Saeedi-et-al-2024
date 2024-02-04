%  V1_cells = Cell_area == 1;
% non_overlapping_V1 = non_overlapping(V1_cells);
% 
% significant_response_any = evoked_any(V1_cells); 
% significant_response_both = both_evoked(V1_cells);
% significant_response_LDG = P_values(V1_cells,2) ;
% significant_response_NCS = P_values(V1_cells,1) ;
% significant_response_DBC = P_values(V1_cells,3) ;
% 
% NCS_PSTH   = grand_psth_ldg1(V1_cells,:);
% LDG_PSTH   = grand_psth_ncs1(V1_cells,:);
% DBC_PSTH   = grand_psth_ctrl1(V1_cells,:);
% PSTH_time_ax = xp(1:end-10);
% 
% NCS_FR = ave_act(V1_cells,1);
% LDG_FR = ave_act(V1_cells,2);
% DBC_FR = ave_act(V1_cells,3);
% mouse_ID = mouse_day_idx(V1_cells,1);
% recording_day = mouse_day_idx(V1_cells,2);
% F1_dominants =logical(F1_dominant) & both_evoked ;
% F1_dominants = F1_dominants(V1_cells);
% 
% Phase_diff=twopi_phase_diff(V1_cells);
% 
% pref_ang_NCS=preferred_ang(V1_cells ,1);
% pref_ang_LDG = preferred_ang(V1_cells ,2);



%% generating Fig 3
load('D:\Data\NCS_2024_figure3.mat')



%% PSTH
figure;
axx=[];alph=.1;
    hold on
    
    idx    =   significant_response_any  & non_overlapping_V1;%& f1f0(:,2)<0
    
    

        mean_ldg=mean(LDG_PSTH(idx,:));
        SE_ldg=std(LDG_PSTH(idx,:))/sqrt(sum(idx));
        h1=errorbar(PSTH_time_ax,mean_ldg,SE_ldg,'r','CapSize',0,'LineWidth',1);
        set(h1.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [h1.Line.ColorData(1:3); 255*alph])
        
        mean_ncs=mean(NCS_PSTH(idx,:));
        SE_ncs=std(NCS_PSTH(idx,:))/sqrt(sum(idx));
        h2=errorbar(PSTH_time_ax,mean_ncs,SE_ncs,'b','CapSize',0,'LineWidth',1);
        set(h2.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [h2.Line.ColorData(1:3); 255*alph])
        
        mean_ctrl=mean(DBC_PSTH(idx,:));
        SE_ctrl=std(DBC_PSTH(idx,:))/sqrt(sum(idx));
        h4=errorbar(PSTH_time_ax,mean_ctrl,SE_ctrl,'Color',[0.4660 0.6740 0.1880],'CapSize',0,'LineWidth',1);
        set(h4.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [h4.Line.ColorData(1:3); 255*alph])

    legend(['LDG, #',num2str(sum(idx))],['NCS, #',num2str(sum(idx))],['DBC, #',num2str(sum(idx))],'Location','northeast')
    legend('boxoff')
    xl=xline(0,'LineWidth',1,'Color','g');
    set(get(get(xl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    xl=xline(1,'LineWidth',1,'Color','g');
    set(get(get(xl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    xlabel('time (s)')
    ax= gca; ax.FontSize = 10; axx=[axx ax];
    ylim([-.05 .5])






%% histogram for angle diff
figure;
idx = significant_response_both  & non_overlapping_V1;
pref_chang = rad2deg(wrapToPi(deg2rad(pref_ang_NCS(idx)-pref_ang_LDG(idx))));
histogram(pref_chang,-186:45:185,'Normalization','probability')
title(['Direction selective cells; N=',num2str(sum(idx))])
% xlabel("LDG_{preferred angle}-NIG_{preferred angle}")
xlabel("changes in preferred angle")
ylabel("probability")
ax=gca;
ax.FontSize = 11;
ax.FontName="Arial";
set(ax,'linewidth',1.5)
box off
axis square
% ylim([0,450])
xticks(-180:90:180)

sum(significant_response_any(non_overlapping_V1))

%% pie chart
sum(significant_response_any(non_overlapping_V1))
V1_ncs_only = sum(significant_response_NCS(non_overlapping_V1) & ~significant_response_LDG(non_overlapping_V1))/sum(significant_response_any(non_overlapping_V1))
V1_ncs_ldg = sum(significant_response_NCS(non_overlapping_V1) & significant_response_LDG(non_overlapping_V1))/sum(significant_response_any(non_overlapping_V1))
V1_ldg_only = sum(~significant_response_NCS(non_overlapping_V1) & significant_response_LDG(non_overlapping_V1))/sum(significant_response_any(non_overlapping_V1))

figure;
pie([V1_ncs_only,V1_ncs_ldg,V1_ldg_only],'%.2f%%')
legend({'ncs only','ncs and ldg','ldg only'})


%% polar histogram for phase difference

phase_diff=deg2rad(Phase_diff(F1_dominants & non_overlapping_V1));
% circ_R= circ_r(b);
% [pval, z] = circ_rtest(b)
[mu, ul, ll] = circ_mean(phase_diff);
% [s, s0] = circ_std(b);

figure
polhist=polarhistogram(phase_diff,25,'EdgeColor','blue','LineWidth',0.2);
hold on
polarplot([0,mu],[0,max(polhist.Values)],'r','LineWidth',3)
title([num2str(sum(F1_dominants)),' units'])

ax = gca;
ax.FontSize = 14; 
ax.FontName="Arial";
set(ax,'linewidth',1)
