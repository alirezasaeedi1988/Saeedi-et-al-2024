


clear
load('D:\Data\Data_figure5.mat')

area_names  = {'V1';'LM';'LI';'LL'};
%% V1 comparison
Ncells  =  cell(1,1);
% figure('units','normalized','outerposition',[0 0 .8 1])
figure;
for area=1%:4
    mean_fr=zeros(2);
    err=zeros(2);
    ci_mean=zeros(2);

    idx =  Cell_area==area &  significant_response_any &  non_overlapping;
    Ncells{area}    = [area_names{area},': ',num2str(sum(idx)),' cells'];
    for i=1:4
        mean_fr(i) = mean(  Firing_rate(idx,i));
        err(i)  = std(  Firing_rate(idx,(i)))/sum(idx);
        pd= fitdist(  Firing_rate(idx,(i)),'Normal');
        ci = pd.paramci-pd.mean;
        ci_mean(i) =ci(2,1);
    end

    % subplot(2,2,area)
    bar(mean_fr)
    hold on
    ngroups = size(mean_fr, 1);
    nbars = size(mean_fr, 2);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, mean_fr(:,i), ci_mean(:,i),'LineWidth', 2,'CapSize',20,'LineStyle', 'none');

    end
    xticklabels({'NCS', 'LDG'})
    legend({'Light off','Light on'},'Location','southeastoutside')
    ylabel({'Firing rate'})
    title(Ncells{area})
    ax = gca;
    ax.FontSize = 14;
    ax.FontName="Arial";
end

%% non intersecting psth
alph = .2;
evoked_type = 'any';
% figure('units','normalized','outerposition',[0 0 0.5 0.6])
figure
hh=[];
for area=1 %:3
    %     subplot(2,2,area);
    hold on
    idx =  Cell_area==area &  significant_response_any &  non_overlapping;

    if sum(idx)>1
        mean_ncs = mean(grand_psth_ncs1(idx,:));
        SE_ncs   = std(grand_psth_ncs1(idx,:))/sqrt(sum(idx));
        h2 = errorbar(xp,mean_ncs,SE_ncs,'b','CapSize',0,'LineWidth',1);
        set(h2.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [h2.Line.ColorData(1:3); 255*alph])

        mean_ldg = mean(grand_psth_ldg1(idx,:));
        SE_ldg   = std(grand_psth_ldg1(idx,:))/sqrt(sum(idx));
        h1 = errorbar(xp,mean_ldg,SE_ldg,'r','CapSize',0,'LineWidth',1);
        set(h1.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [h1.Line.ColorData(1:3); 255*alph])


        mean_ncs = mean(opto_psth_ncs(idx,:));
        SE_ncs   = std(opto_psth_ncs(idx,:))/sqrt(sum(idx));
        h2 = errorbar(xp,mean_ncs,SE_ncs,'c','CapSize',0,'LineWidth',1);
        set(h2.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [h2.Line.ColorData(1:3); 255*alph])

        mean_ldg = mean(opto_psth_ldg(idx,:));
        SE_ldg   = std(opto_psth_ldg(idx,:))/sqrt(sum(idx));
        h1 = errorbar(xp,mean_ldg,SE_ldg,'m','CapSize',0,'LineWidth',1);
        set(h1.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [h1.Line.ColorData(1:3); 255*alph])


        mean_ctrl = mean(grand_psth_ctrl1(idx,:));
        mean_ctrl = mean_ctrl-mean(mean_ctrl);
        SE_ctrl   = std(grand_psth_ctrl1(idx,:))/sqrt(sum(idx));
        h2 = errorbar(xp,mean_ctrl,SE_ctrl,'k','CapSize',0,'LineWidth',1);
        set(h2.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [h2.Line.ColorData(1:3); 255*alph])

    end

    %     legend(['NCS, #',num2str(sum(idx))],['LDG, #',num2str(sum(idx))],...
    %         ['NCS+L, #',num2str(sum(idx))],['LDG+L, #',num2str(sum(idx))],'Location','northeast')
    if area==1
        legend('NCS','LDG','NCS+L','LDG+L','DBC','Location','northeast')
        legend('boxoff')
    end


    xl = xline(0,'LineWidth',1,'Color','g');
    set(get(get(xl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    xl = xline(1,'LineWidth',1,'Color','g');
    set(get(get(xl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    xlabel('time (s)')
    title([area_names{area},' #',num2str(sum(idx))])
    ax = gca; ax.FontSize = 10;ax.FontName = "Arial";
    %     hh = [hh ax];
    %     ylim([-.5,.8])
end
