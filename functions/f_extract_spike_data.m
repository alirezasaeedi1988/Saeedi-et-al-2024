function f_extract_spike_data(myKsDir)
%%% extracting spike times and...
%%% loading files

inputfile_stm = [myKsDir,myKsDir(end-15:end),'_stimulus_time.mat'];
load(inputfile_stm)
spike_times   = readNPY([myKsDir ,'\spike_times.npy']);
spike_cluster = readNPY([myKsDir ,'\spike_clusters.npy']);
[cids, cgs]   = readClusterGroupsCSV([myKsDir ,'\cluster_group.tsv']);

if length(stimulus_point)==3 %#ok<*USENS>
    session_numbers=2:3;
    session_names={'RF';'size';'neon'};
elseif length(stimulus_point)==4 
    session_numbers=2:3;
    session_names={'RF';'size';'neon';'MAE'};
elseif length(stimulus_point)==5
    session_numbers=2:4;
    session_names={'RF';'size';'neon';'size2';'MAE'};
elseif length(stimulus_point)==6
    session_numbers=2:5;
    session_names={'RF';'size';'contextual';'size2';'neon';'MAE'};
end

 for session_number=session_numbers
    
    session_name=session_names{session_number};

    fs  =30000;                                          %sampling frequency
    trial_schedule=stimulus_schedule{session_number};
    tf  =stimulus_TemporalFreq{session_number};          % temporal frequency (cycle/s)
    Nt=min([length(trial_schedule),length(stimulus_time{session_number})-2]);          %#ok<*IDISVAR>
    pst=0.2;                                             % pre-stimulus time sec
    dx= 1e-3;                                            % bin size
    stim_len=mean(stimulus_length{session_number})/fs;
%     band_test=-pst:pst:stim_len;
    %%%
    Spt=floor(stimulus_point{session_number}-(pst*fs));      %starting point of trials in the data
    nsample=floor(mean(diff(Spt(2:end-1))));                 %number of samples for each trials
    Sps=stimulus_point{session_number};                      %starting point of stimulus
    %%%
    if strcmp(session_name,'neon')
        pst=0.3;                                                 % pre-stimulus time sec
        Spt=floor(stimulus_point{session_number}-(pst*fs));      %starting point of trials in the data
        nsample=floor(mean(diff(Spt)));                          %number of samples for each trials
        trial_schedule(Nt+1:end,:)=[];
        if min(size(trial_schedule))==2
            realcolor_loc=trial_schedule(:,1)==2;
            trial_schedule(realcolor_loc,2)=trial_schedule(realcolor_loc,2).*10;
            trial_schedule=trial_schedule(:,2);
        else
            phys1_loc=(trial_schedule(:,1)==1 & trial_schedule(:,2)==2);
            phys2_loc=(trial_schedule(:,1)==2 & trial_schedule(:,2)==2);
            ctrl_loc=(trial_schedule(:,1)==2 & trial_schedule(:,2)==1);
            
            trial_schedule(phys1_loc,3)=trial_schedule(phys1_loc,3).*10;
            trial_schedule(phys2_loc,3)=trial_schedule(phys2_loc,3).*100;
            trial_schedule(ctrl_loc,3)=trial_schedule(ctrl_loc,3).*1000;
            
            trial_schedule=trial_schedule(:,3);
        end
    elseif strcmp(session_name,'contextual')
        number_of_discarding_trials=length(trial_schedule)-length(Sps);
        trial_schedule(1:number_of_discarding_trials)=[]; 
    end
    band_test=[-pst 0 stim_len/2];
    condLabels = unique(trial_schedule);
    %%%
    xp=-pst:dx:(nsample/fs)-pst;
    [~, onset_ind] = min(abs(xp));
    [~, resp_ind]  = min(abs(xp-0.05));
    [~,offset_ind] = min(abs(xp-stim_len));
%     [~,onset_resp_ind]=min(abs(xp-(pst/2)));
    %%%
    good_cid= cgs==2;
    good_cid=cids(good_cid);
    spikes=cell(1,2);
    pdfs=cell(1,2);
    RFs=cell(1,2);    
    smoothed_pdf=cell(1,2);
    evoked_cids=cell(1,1);
    average_latency=cell(1,2);
    trial_per_cond=cell(1,length(good_cid));
    e=1;
    ep=1;
    for neuron=1:length(good_cid) %%%loop over clusters, spike time for selected cluster is extracted
        RF=zeros(length(condLabels),length(xp));
        hist0=zeros(length(condLabels),length(xp));
        smoothed_hist=zeros(length(condLabels),length(xp));
        cid= spike_cluster==good_cid(neuron);
        Cluster=spike_times(cid);
        Cluster=Cluster';
        Cluster = double(Cluster);
        con_Ntrial=zeros(1,length(condLabels));
        test_vector=[];
        spike=[];
        for cond=1:length(condLabels)%%% loop over all condition for selected neuron
            ind_trGroup= trial_schedule(1:Nt)==condLabels(cond);
            cond_spt=Spt(ind_trGroup);
            cond_sps=Sps(ind_trGroup);
            con_Ntrial(cond)=length(cond_spt);
            Trial=cell(con_Ntrial(cond),1);
            
            for trial_number=1:con_Ntrial(cond) %% loop over all trails in a condition for alignment spike times
                ind_tr= cond_spt(trial_number)<=Cluster & Cluster<cond_spt(trial_number)+nsample;
                Tri=Cluster(ind_tr);
                if ~isempty(Tri)
                    %%% align
                    Tri=(Tri-cond_sps(trial_number))./(fs);
%                     hhh=histogram(Tri,band_test);
%                     test_vector=[test_vector;hhh.Values];
                    spike=[spike;Tri' , ones(length(Tri),1)*trial_number , ones(length(Tri),1)*cond]; %#ok<AGROW>
                    Trial{trial_number}=Tri;
                    hhh=zeros(1,length(band_test)-1);
                    for iii=1:length(hhh)
                        hhh(iii)=sum(Tri>=band_test(iii)& Tri<band_test(iii+1))/(band_test(iii+1)-band_test(iii));
                    end
                    test_vector=[test_vector;hhh]; %#ok<AGROW>

                end
            end
            %%%PSTH calculation
            P=zeros(1,length(xp));
            for trial_number=1:con_Ntrial(cond) %% loop over all trails in a condition for psth
                if ~isempty(Trial{trial_number})
                    for jj=1:length(Trial{trial_number})
                        n=floor(((Trial{trial_number}(jj)-xp(1))/dx))+1;
                        P(n)=P(n)+1;
                    end
                end
            end
            PP=P./(con_Ntrial(cond)*dx); %% normalizing psth with number of trials and bin size
            RF(cond,:)=P;        %relative frequency
            hist0(cond,:)=PP;
            smoothed_hist(cond,:)=smoothdata(PP,'gaussian',15);%%% In 'smoothdata' the SD is fixed to 1/5th of window width
            
        end
        %%% poisson test for evoked and non-evoked neuron
        y                =  sum(RF);
        teta=f_resp_latency(y,onset_ind,offset_ind);
        if ~isempty(teta)
            p_evoked_cids{ep}=good_cid(neuron);
            average_latency{ep}=teta;
            ep=ep+1;
        end
        %%% wilcoson test for evoked and non-evoked neuron
        if ~isempty(test_vector)
            p_value=zeros(1,length(hhh)-1);
            h_value=zeros(1,length(hhh)-1);
            for iii=1:length(hhh)-1
                [p_value(iii),h_value(iii)]=ranksum(test_vector(:,1),test_vector(:,iii+1),'alpha', 0.001);
            end
            if h_value(1)==1
                spikes{e}=spike;
                trial_per_cond{e}=con_Ntrial;
                pdfs{e}=hist0;
                RFs{e}=RF;
                smoothed_pdf{e}=smoothed_hist;
                evoked_cids{e}=good_cid(neuron);
                e=e+1;
            end
        end
    end
    if strcmp(session_name,'neon')
        sr=f_surround_suppression(myKsDir);
        spike_data.surround_suppression=sr;
    end
    spike_data.spikes          = spikes;
    spike_data.RFs             = RFs;
    spike_data.pdfs            = pdfs;
    spike_data.smoothed_pdf    = smoothed_pdf;
    spike_data.pst             = pst;
    spike_data.dx              = dx;
    spike_data.xp              = xp;
    spike_data.cluster_id      = good_cid;
    spike_data.trial_per_cond  = trial_per_cond;
    spike_data.stim_len        = stim_len;
    spike_data.evoked_cids     = evoked_cids;
    spike_data.average_latency = average_latency;
    spike_data.p_evoked_cids   = p_evoked_cids;
    save([myKsDir,'\',session_name,'_spike_data.mat'],'spike_data')
    
end