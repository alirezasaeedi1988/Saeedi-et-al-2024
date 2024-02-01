function f_neon_tuning_curve(myKsDir)

session_name='neon';

outputDir           = fullfile(myKsDir,session_name);
    if ~exist(outputDir, 'dir')
       mkdir(outputDir)
    end

load([myKsDir,'\',session_name,'_spike_data']) %#ok<LOAD>
pdfs=spike_data.pdfs;
xp=spike_data.xp;
% cid=spike_data.cluster_id;
cid=spike_data.evoked_cids;

stim_len=spike_data.stim_len;



tuningYlabel='Spike/sec';
% tuningYlabel='Average firing rate';
tuningXlabel='Grating angles';

stim_angle=[0 45 90 135 180 225 270 315];
stim_type=[1 2];
neon_type=[1 2];
Ncond=length(stim_type)*length(stim_angle)*length(neon_type);     %number of condition

heatmapDim1=length(stim_angle);
heatmapDim2=length(stim_type)*length(neon_type);
       

time_window=[0 stim_len/2];                 % averaging time  window 




win_ind = find(xp>=time_window(1) & xp<time_window(2));

Nc=length(pdfs);       %number of clusters

for neuron=1:Nc %%%loop over clusters
    heatmap=zeros(heatmapDim1,heatmapDim2);
    histo=pdfs{neuron};
    for cond=1:Ncond%%% loop over all condition for selected cluster
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% average firing rate
        heatmap(cond)=mean(histo(cond,win_ind(1):win_ind(end)));
    end
    heatmap1=zeros(heatmapDim1,heatmapDim2);
    heatmap1(:,1)=smoothdata(heatmap(:,1),'gaussian');%%% In 'smoothdata' the SD is fixed to 1/5th of window width
    heatmap1(:,2)=smoothdata(heatmap(:,2),'gaussian');%%% In 'smoothdata' the SD is fixed to 1/5th of window width
    heatmap1(:,4)=smoothdata(heatmap(:,4),'gaussian');%%% In 'smoothdata' the SD is fixed to 1/5th of window width
    
    figure('units','normalized','outerposition',[0 0 0.5 0.6]);set(gcf, 'Visible', 'off');
    plot(heatmap(:,1),'bd','LineWidth',2);hold on;
    plot(heatmap1(:,1),'-.b','LineWidth',2);
    plot(heatmap(:,2),'r*','LineWidth',2);
    plot(heatmap1(:,2),'-r','LineWidth',2);
    
%     plot(heatmap(:,3),'-m*','LineWidth',2);
    plot(heatmap(:,4),'kd','LineWidth',2);
    plot(heatmap1(:,4),'-k','LineWidth',2);
    
    box off
%     legend('neon-color','physical-1','physical-2','ctrl','Location','northwest')
    legend('NIG','smoothed NIG','LDG','smoothed LDG','ctrl','smoothed ctrl','Location','northwest')
    legend('boxoff')
    title(['direction tuning curve cluster:',num2str(cid{neuron})]);
    xlabel(tuningXlabel);ylabel(tuningYlabel);
    %     set(gca,'XTick',stim_angle)
    set(gca,'xticklabel',stim_angle)
    set(gca,'linewidth',1.5)
    ax = gca;
    ax.FontSize = 18;
    ax.FontName = "Arial";
    saveas(gcf,[outputDir,'\',session_name,'_tuning curve_cluster-',num2str(cid{neuron}),'.png'])
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperOrientation','landscape');
    saveas(gcf,[outputDir,'\',session_name,'_tuning curve_cluster-',num2str(cid{neuron}),'.pdf'])
    close gcf

end
    
