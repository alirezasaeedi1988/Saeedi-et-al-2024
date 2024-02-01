function f_neon_conditional_raster_psth_(myKsDir)

session_name='neon';

outputDir           = fullfile(myKsDir,session_name);
    if ~exist(outputDir, 'dir')
       mkdir(outputDir)
    end

load([myKsDir,'\',session_name,'_spike_data.mat'])
smoothed_pdf=spike_data.smoothed_pdf;
xp=spike_data.xp;
cid=spike_data.evoked_cids;
spikes=spike_data.spikes;
stim_len=spike_data.stim_len;
trial_per_cond=spike_data.trial_per_cond;


%% plot options

        Nhp=8; Nvp=4;
        raster_loc=[1 3 9 11 17 19 25 27];raster_loc=[raster_loc raster_loc+1];
        psth_loc=[5 7 13 15 21 23 29 31];psth_loc=[psth_loc psth_loc];
        tuningYlabel='Average firing rate';
        tuningXlabel='Grating angles';

stim_angle=[0 45 90 135 180 225 270 315];
stim_type=[1 2];
Ncond=length(stim_type)*length(stim_angle);                 %number of condition
heatmapDim1=length(stim_angle);
heatmapDim2=length(stim_type);    
time_window=[0 0.6];                 % averaging time  window 
win_ind = find(xp>=time_window(1) & xp<time_window(2));
Nc=length(smoothed_pdf);       %number of clusters

for neuron=1:Nc %%%loop over clusters
    heatmap=zeros(heatmapDim1,heatmapDim2);
    spike_in_neuron=spikes{neuron};
    if ~isempty(spike_in_neuron)
        figure('units','normalized','outerposition',[0 0 1 1]);
        set(gcf, 'Visible', 'off');
        histo=smoothed_pdf{neuron};
        for cond=1:Ncond%%% loop over all condition for selected cluster
            
            spike_in_cond=spike_in_neuron(spike_in_neuron(:,3)==cond,1:2);
            subplot(Nhp,Nvp,raster_loc(cond))
            for k=1:trial_per_cond{neuron}(cond) %% loop over all trails in a condition for raster plot
                Tri=spike_in_cond(spike_in_cond(:,2)==k,1);
                if ~isempty(Tri)
                    y=ones(1,length(Tri)).*(k);
                    scatter(Tri,y,1,'black','fill')
                    hold on
                end
            end
            set(gca,'xticklabel',[]);
            set(gca,'YTick',[]);
            xlim([xp(1) xp(end)]);
            ylim([0 trial_per_cond{neuron}(cond)])
            if cond<9
                ylabel('neon color')
            else
                ylabel('physical')
            end
            %%%psth plot
            subplot(Nhp,Nvp,psth_loc(cond):psth_loc(cond)+1)
            hold on
            sp=smooth(histo(cond,:),'lowess');
%             sp=histo(cond,:);
            if cond<9
                plot(xp,sp,'-.b')  
                ylabel(['angle=',num2str(stim_angle(cond))])
                
            else
                plot(xp,sp,'-r');
                legend('neon color grating','physical grating')
                rectangle('Position',[0 0 stim_len max(histo(:))],'EdgeColor','g')

            end
%             set(gca,'xticklabel',[]);
            set(gca,'YTick',[]);
            xlim([xp(1) xp(end)]);
            ylim([min(histo(:)) max(histo(:))])
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% average firing rate
            
            heatmap(cond)=mean(histo(cond,win_ind(1):win_ind(end)));
        end
        saveas(gcf,[outputDir,'\',session_name,'_conditional PSTH_cluster-',num2str(cid{neuron}),'.png'])
        close gcf
    end
    figure('units','normalized','outerposition',[0 0 0.5 0.6]);set(gcf, 'Visible', 'off');
    plot(heatmap(:,1),'-.bd');
    hold on
    plot(heatmap(:,2),'-r*');
    
    legend('neon color grating','physical grating','Location','northwest')
    title(['neon color spread ',num2str(cid{neuron}),', window=',num2str(time_window(1)*1000),'-',num2str(time_window(2)*1000),'ms']);
    xlabel(tuningXlabel);ylabel(tuningYlabel);
%     set(gca,'XTick',stim_angle)
        set(gca,'xticklabel',stim_angle)
    saveas(gcf,[outputDir,'\',session_name,'_tuning curve_win-',num2str((time_window(2)-time_window(1))*1000),'ms_cluster-',num2str(cid{neuron}),'.png'])
    close gcf

end
    
