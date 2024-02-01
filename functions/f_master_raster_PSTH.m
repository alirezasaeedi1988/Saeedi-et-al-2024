function f_master_raster_PSTH(myKsDir)


% master raster plot
% myKsDir=KS_dirs{day}

addpath(genpath('D:\matwork\KiloSort'))
addpath(genpath('D:\matwork\npy-matlab'))

inputfile_stm           = [myKsDir,myKsDir(end-15:end),'_stimulus_time.mat'];
load(inputfile_stm)
spike_times=readNPY([myKsDir ,'\spike_times.npy']);
spike_cluster=readNPY([myKsDir ,'\spike_clusters.npy']);




fs=30000; %sampling frequency
pst=0.3;                                   %pre-stimulus time sec
dx= 20e-3;


if length(stimulus_point)==3  %#ok<*USENS>
    session_numbers=2:3;
    session_names={'RF';'size';'neon'};
elseif length(stimulus_point)==4 
    session_numbers=2:3;
    session_names={'RF';'size';'neon';'MAE'};
elseif length(stimulus_point)==5
    session_numbers=2:4;
    session_names={'RF';'size';'neon';'size2';'MAE'};
%     session_names={'RF';'size';'contextual';'size2';'neon'};
elseif length(stimulus_point)==6
    session_numbers=2:5;
    session_names={'RF';'size';'contextual';'size2';'neon';'MAE'};
end





for session_number=session_numbers %%%% loop over sessions
    
    session_name=session_names{session_number};
    load([myKsDir,'\',session_name,'_spike_data.mat'])
    evoked_cid=spike_data.evoked_cids;
    
    if strcmp(session_name,'neon')
        pst=0.3;
    end
    outputDir           = fullfile(myKsDir,session_name);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir)
    end

    Nt=length(stimulus_point{session_number});  %#ok<*IDISVAR>
    spt=floor(stimulus_point{session_number}-(pst*fs));         %starting point of trials in the data
    nsample=floor(mean(diff(spt(2:end-2))));            %number of samples for each trials
%     time=(0:nsample)./fs-pst;

for neuron=1:length(evoked_cid)%%%loop over clusters
    
    cid= spike_cluster==evoked_cid{neuron};
    spikes_in_cluster=spike_times(cid);
    spikes_in_cluster=spikes_in_cluster';
    spikes_in_cluster = double(spikes_in_cluster);
    Max=zeros(1,Nt-1);
    Min=zeros(1,Nt-1);
    Trial=cell(Nt,1);
    figure;set(gcf, 'Visible', 'off');
    subplot(2,1,1)
    for k=1:Nt-1
        ind_tr= spt(k)<=spikes_in_cluster & spikes_in_cluster<=spt(k)+nsample;
        Tri=spikes_in_cluster(ind_tr);
        %%% align
        if ~isempty(Tri)
            Tri=(Tri-stimulus_point{session_number}(k))./(fs);
            y=ones(1,length(Tri)).*(k);
            scatter(Tri,y,0.15,'black','fill')
            hold on
            Trial{k}=Tri;
            Max(k)= max(Tri);
            Min(k)= min(Tri);
        end
        
    end
    rectangle('Position',[0 0 mean(stimulus_length{session_number})/fs Nt],'EdgeColor','r')
%     set(gca,'XTick',[]);
%     xlabel('time(second)')
    ylabel('trials')
    title(['cluster id: ',num2str(evoked_cid{neuron})]);
    ylim([0 Nt])
    set(gca,'XTick',-pst:0.1:max(Max));set(gca,'yTick',0:floor(Nt/5):Nt)
%%%PSTH plot
   
    xp=min(Min):dx:max(Max);
%     xp=[xp,xp(end)+1];
    p_len=length(xp);%floor((max(Max)-min(Min))/dx)+1;
    P=zeros(1,p_len);
    for k=1:Nt
        if ~isempty(Trial{k})
            for jj=1:length(Trial{k})
                n=floor(((Trial{k}(jj)-min((Min)))/dx))+1;
                P(n)=P(n)+1;
            end
        end
        
    end
    P=P./(Nt*dx); %% normalizing psth with number of trials and bin size
    [~, zero_ind]=min(abs(xp));
    pre_stimulus_baseline=mean(P(1:zero_ind-1));
    P=(P-pre_stimulus_baseline)/pre_stimulus_baseline; %% normalizing psth by pre_stimulus_baseline
    subplot(2,1,2)
    plot(xp,P)
    xlabel('time(second)')
    ylabel('normalized firing rate')
    legend('PSTH')
    set(gca,'XTick',-pst:0.1:max(Max))
%     savefig(['rasterPlot',num2str(i),'.fig'])
    saveas(gcf,[outputDir,'\',session_name,'_raster_PSTH_cluster_',num2str(evoked_cid{neuron}),'.png'])
    close gcf
end
end   
