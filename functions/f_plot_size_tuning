function f_plot_size_tuning(myKsDir)
% specific size tunning


% myKsDir                 = 'D:\matwork\Data\M145_2019_07_12';

session_name='size';

outputDir           = fullfile(myKsDir,session_name);
    if ~exist(outputDir, 'dir')
       mkdir(outputDir)
    end
    
load([myKsDir,'\',session_name,'_spike_data.mat'])
inputfile_stm = [myKsDir,myKsDir(end-15:end),'_stimulus_time.mat'];
load(inputfile_stm)

pdfs=spike_data.pdfs;
xp=spike_data.xp;
cid=spike_data.evoked_cids;

% stim_size=[5:5:30 35:10:55];
% stim_size=[4 6 8 10:5:30 35:10:55 60];
% stim_sizes= [2.5 5:5:45];  %[2.5 5:5:45 55 65];
stim_sizes = stimulus_sizes{2};
stim_angle=stimulus_angles{2};
Ncond=length(stim_angle)*length(stim_sizes);                 %number of condition

%% plot options

tuningYlabel='Spike/sec';
tuningXlabel='Stimulus size(deg)';
heatmapDim1=length(stim_sizes);
heatmapDim2=length(stim_angle);
       

time_window=[0 0.4];                        % averaging time  window 

%%%
win_ind = find(xp>=time_window(1) & xp<time_window(2));

Nc=length(pdfs);       %number of clusters


for neuron=1:Nc
    heatmap=zeros(heatmapDim1,heatmapDim2);
    histo=pdfs{neuron};
   
    for cond=1:Ncond%%% loop over all condition for selected cluster
        
        heatmap(cond)=mean(histo(cond,win_ind(1):win_ind(end)));
            
    end
    f1 = fit(stim_sizes',heatmap(:,1),'smoothingspline');
    f2 = fit(stim_sizes',heatmap(:,2),'smoothingspline');
%     heatmap(:,1)=smooth(heatmap(:,1));
%     heatmap(:,2)=smooth(heatmap(:,2));
    figure('units','normalized','outerposition',[0 0 0.5 0.6]);set(gcf, 'Visible', 'off');
    hL1=plot(f1,stim_sizes,heatmap(:,1),'bd');
    hold on
    box off
    set(get(get(hL1(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hL2=plot(f2,stim_sizes,heatmap(:,2),'r*');
    set(get(get(hL2(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(hL1,'LineWidth',1)
    set(hL2,'LineWidth',1)
    set(hL1(2),'color','b')
    set(hL1,'MarkerFaceColor','b','MarkerSize',5)
    legend('Vertical grating','Horizontal grating','Location','northeast')
    legend('boxoff')
    title(['Size tuning for neuron ',num2str(cid{neuron}),', window=',num2str(time_window(1)*1000),'-',num2str(time_window(2)*1000),'ms']);
    xlabel(tuningXlabel);ylabel(tuningYlabel);
    set(gca,'XTick',stim_sizes)
    set(gca,'linewidth',1.5)
    ax=gca;
    ax.FontSize = 18;
    ax.FontName = "Arial";
%     set(gca,'xticklabel',stim_size)
    saveas(gcf,[outputDir,'\',session_name,'_win-',num2str((time_window(2)-time_window(1))*1000),'ms_cluster-',num2str(cid{neuron}),'.png'])
    set(gcf,'PaperPositionMode','auto'); 
    set(gcf,'PaperOrientation','landscape');
    saveas(gcf,[outputDir,'\',session_name,'_win-',num2str((time_window(2)-time_window(1))*1000),'ms_cluster-',num2str(cid{neuron}),'.pdf'])
    close gcf
 
end
    
