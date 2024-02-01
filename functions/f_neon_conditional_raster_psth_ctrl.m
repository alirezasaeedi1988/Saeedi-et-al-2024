function f_neon_conditional_raster_psth_ctrl(myKsDir)

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

        Nhp=12; Nvp=4;
        raster_loc=[1 3 13 15 25 27 37 39];raster_loc=[raster_loc raster_loc+1 raster_loc+5 raster_loc+4];
        psth_loc=[9 11 21 23 33 35 45 47];psth_loc=[psth_loc psth_loc psth_loc psth_loc];
        

stim_angle=[0 45 90 135 180 225 270 315];
stim_type=[1 2];
neon_type=[1 2];
Ncond=length(stim_type)*length(stim_angle)*length(neon_type);                 %number of condition


Nc=length(smoothed_pdf);       %number of clusters

for neuron=1:Nc %%%loop over clusters
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
                ylabel('neon')
            elseif (8<cond)&& (cond<17)
                ylabel('phy-1')
                
            elseif (16<cond)&& (cond<25)
                ylabel('phy-2')
            elseif cond>24
                ylabel('ctrl')
            end
            %%%psth plot
            subplot(Nhp,Nvp,psth_loc(cond):psth_loc(cond)+1)
            hold on
            sp=histo(cond,:);
            if cond<9
                plot(xp,sp,'-.b')
                ylabel(['angle=',num2str(stim_angle(cond))])
            elseif (8<cond)&& (cond<17)
                plot(xp,sp,'-r')
            elseif (16<cond)&& (cond<25)
                plot(xp,sp,'-m')
            elseif cond>24
                plot(xp,sp,'-k');
                legend('neon-color','physical-1','physical-2','ctrl')
                rectangle('Position',[0 0 stim_len max(histo(:))],'EdgeColor','g')

            end
            set(gca,'xticklabel',[]);
            set(gca,'YTick',[]);
            xlim([xp(1) xp(end)]);
            ylim([min(histo(:)) max(histo(:))])
                        
        end
        saveas(gcf,[outputDir,'\',session_name,'_conditional PSTH_cluster-',num2str(cid{neuron}),'.png'])
        close gcf
    end
end
    
