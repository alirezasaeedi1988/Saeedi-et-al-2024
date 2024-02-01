function f_neon6_plot_pref_raster_smooth(myKsDir,cond_preferred_angle,psth_locs,spike_data,neuron,...
    neuron_latency5p,neuron_latency3,bin_size,area_name,session_name,f1f0,f2f0,p_values)

neuron_latency5p=neuron_latency5p./1000;
neuron_latency3=neuron_latency3./1000;


% session_name='neon';

outputDir           = fullfile(myKsDir,session_name);
if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end

% load([myKsDir,'\',session_name,'_spike_data.mat'])
smoothed_pdf=spike_data.smoothed_pdf{neuron};

xp=spike_data.xp;
spike_in_neuron=spike_data.spikes{neuron};
stim_len=spike_data.stim_len;
Ntrial_per_cond=spike_data.trial_per_cond{neuron}(psth_locs);
id=spike_data.evoked_cids{neuron};
% Ntrial_per_cond=trial_per_cond{neuron};
% spike_in_neuron= spikes{neuron};
% smoothed_pdf=smoothed_pdfs{neuron};
%% plot options

Nhp=3; Nvp=2;
% raster_loc=[1 2 6 5 3 4];




if ~isempty(spike_in_neuron)
    figure('units','normalized','outerposition',[0 0 0.4 .75],'Visible','off');
    histo=smoothed_pdf(psth_locs,:);
    MAX1=[];
    for cond=1:4%%% loop over all condition for selected cluster and preferred angle
        spike_in_cond=spike_in_neuron(spike_in_neuron(:,3)==psth_locs(cond),1:2);
        subplot(Nhp,Nvp,cond)
        scatter(spike_in_cond(:,1),spike_in_cond(:,2),2,'black','fill') %try this instead of loop
        
        %         set(gca,'xticklabel',[]);
        %         set(gca,'YTick',[]);
        xlim([xp(1) xp(end)]);
        ylim([0 Ntrial_per_cond(cond)])
        ylabel('trial')
        xline(0,'g','LineWidth',1);
        xline(1,'g','LineWidth',1);
%         rectangle('Position',[0 0 stim_len Ntrial_per_cond(cond)],'EdgeColor','g','LineWidth',2)
        if ~isnan(neuron_latency5p(cond))
            xline(neuron_latency5p(cond),'-.r',{'L5%',num2str(neuron_latency5p(cond))},'LineWidth',1,'LabelVerticalAlignment','top','LabelOrientation','horizontal');
        end
        if p_values(cond)<0.01 && neuron_latency3(cond)<.3
            xline(neuron_latency3(cond),'-.b',{'L3',num2str(neuron_latency3(cond))},'LineWidth',1,'LabelVerticalAlignment','bottom','LabelOrientation','horizontal');
        end
        switch cond
            case 1
                title({['NCS ','f1f0=',num2str(f1f0(cond),'%4.2f'),' f2f0=',num2str(f2f0(cond),'%4.2f')];[' ang.=',num2str(cond_preferred_angle(cond)),' Pv=',num2str(p_values(cond),'%2.4f')]})
            case 2
                title({['LDG ','f1f0=',num2str(f1f0(cond),'%4.2f'),' f2f0=',num2str(f2f0(cond),'%4.2f')];[' ang.=',num2str(cond_preferred_angle(cond)),' Pv=',num2str(p_values(cond),'%2.4f')]})
            case 3
                title({['NCS+ring ','f1f0=',num2str(f1f0(cond),'%4.2f'),' f2f0=',num2str(f2f0(cond),'%4.2f')];[' ang.=',num2str(cond_preferred_angle(cond)),' Pv=',num2str(p_values(cond),'%2.4f')]})
            case 4
                title({['NCS+square ','f1f0=',num2str(f1f0(cond),'%4.2f'),' f2f0=',num2str(f2f0(cond),'%4.2f')];[' ang.=',num2str(cond_preferred_angle(cond)),' Pv=',num2str(p_values(cond),'%2.4f')]})
            case 5
                title({['NCS+circle ','f1f0=',num2str(f1f0(cond),'%4.2f'),' f2f0=',num2str(f2f0(cond),'%4.2f')];[' ang.=',num2str(cond_preferred_angle(cond)),' Pv=',num2str(p_values(cond),'%2.4f')]})
            case 6
                title({['LDG+circle ','f1f0=',num2str(f1f0(cond),'%4.2f'),' f2f0=',num2str(f2f0(cond),'%4.2f')];[' ang.=',num2str(cond_preferred_angle(cond)),' Pv=',num2str(p_values(cond),'%2.4f')]})
        end
        ax = gca;
        ax.FontSize = 7;
        
        %%%psth plot
        subplot(Nhp,Nvp,5:6)
        hold on
        sp=histo(cond,:);
        s1=length(sp);
        m  = s1 - mod(s1, bin_size);
        sp=mean(reshape(sp(1:m),bin_size,[]));
        xpp=mean(reshape(xp(1:m),bin_size,[]));
        switch cond
            case 1
                plot(xpp,sp,'-b','LineWidth',1)
            case 2
                plot(xpp,sp,'-r','LineWidth',1)
            case 3
                plot(xpp,sp,'-y','LineWidth',1)
            case 4
                plot(xpp,sp,'-k','LineWidth',1);
            case 5
                plot(xpp,sp,'-c','LineWidth',1);
            case 6
                plot(xpp,sp,'-m','LineWidth',1);
        end
        MAX1=[MAX1 max(sp(:))];
    end
    ylabel('Firing rate')
    legend('NCS','LDG','NCS+ring block','NCS+square block')
    rectangle('Position',[0 0 stim_len max(MAX1)],'EdgeColor','g','LineWidth',.75)
    %         xticks(xp(1):.05:xp(end))
    %         set(gca,'xticklabel',[]);
    %         set(gca,'YTick',[]); 
    xlim([xp(1) (xp(end)+xp(end)/1.5)]);
    ylim([min(histo(:)) max(MAX1)])
    xlabel('time (second)')
    sgtitle([area_name,' neuron:',num2str(id)],'FontSize',8)
    ax = gca;
    ax.FontSize = 7;
    saveas(gcf,[outputDir,'\',num2str(bin_size),session_name,'_preferred raster_cluster-',num2str(id),'.pdf'])
    saveas(gcf,[outputDir,'\',num2str(bin_size),session_name,'_preferred raster_cluster-',num2str(id),'.jpg'])

    %% save in different folders
    if ~isempty(session_name)

    if p_values(1)<0.001 && p_values(2)<0.001
    saveas(gcf,[myKsDir(1:end-15),'pv0.001\',myKsDir(end-14:end),'_cluster-',num2str(id),'.pdf'])
    end
    if p_values(1)<0.005 && p_values(2)<0.005
        saveas(gcf,[myKsDir(1:end-15),'pv0.005\',myKsDir(end-14:end),'_cluster-',num2str(id),'.pdf'])
    end
    if p_values(1)<0.01 && p_values(2)<0.01
        saveas(gcf,[myKsDir(1:end-15),'pv0.01\',myKsDir(end-14:end),'_cluster-',num2str(id),'.pdf'])
    end    
    if p_values(1)<0.001 && p_values(2)<0.001 && f1f0(2)<0
    saveas(gcf,[myKsDir(1:end-15),'pv0.001_complex\',myKsDir(end-14:end),'_cluster-',num2str(id),'.pdf'])
    end
    if p_values(1)<0.005 && p_values(2)<0.005 && f1f0(2)<0
        saveas(gcf,[myKsDir(1:end-15),'pv0.005_complex\',myKsDir(end-14:end),'_cluster-',num2str(id),'.pdf'])
    end
    if p_values(1)<0.01 && p_values(2)<0.01 && f1f0(2)<0
        saveas(gcf,[myKsDir(1:end-15),'pv0.01_complex\',myKsDir(end-14:end),'_cluster-',num2str(id),'.pdf'])
    end
    
    if p_values(1)<0.001 && p_values(2)<0.001 && f1f0(2)<0 && f1f0(1)<0
    saveas(gcf,[myKsDir(1:end-15),'pv0.001_both_complex\',myKsDir(end-14:end),'_cluster-',num2str(id),'.pdf'])
    end
    if p_values(1)<0.005 && p_values(2)<0.005 && f1f0(2)<0 && f1f0(1)<0
        saveas(gcf,[myKsDir(1:end-15),'pv0.005_both_complex\',myKsDir(end-14:end),'_cluster-',num2str(id),'.pdf'])
    end 
    if p_values(1)<0.01 && p_values(2)<0.01 && f1f0(2)<0 && f1f0(1)<0
        saveas(gcf,[myKsDir(1:end-15),'pv0.01_both_complex\',myKsDir(end-14:end),'_cluster-',num2str(id),'.pdf'])
    end    
    end
    close gcf
end


