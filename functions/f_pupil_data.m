function f_pupil_data(myKsDir,BlDir,Video_folder)

%% pupil analysis
%% this function analyzes the DeepLabCut results of pupil detection and syncronizes them with stiimulus onsets!!!



% Video_folder = 'Y:\neon_pupil-Alireza-2022-09-01\laser_neon';
anlyized_video_list = dir(fullfile(Video_folder,'*.mp4'));

inputfile_stm = [myKsDir,myKsDir(end-15:end),'_stimulus_time.mat'];
load(inputfile_stm)


if length(stimulus_point)==5%#ok<*USENS>
    session_numbers=[ 5 ];  %% add analysis of opto here too
    session_names = {'RF';'opto';'rand';'RF2';'neon'};
    Nsessions = length(session_numbers);
    %% 
    sessions_strct=dir(fullfile(BlDir));
    sessions_dirindx=zeros(1,Nsessions);
    jj=1;
    for k=1:length(sessions_strct)
        if length(sessions_strct(k).name)==19 && strcmp(sessions_strct(k).name(1:2),'20')
            %keyboard
            sessions_dirindx(jj)= k;
            jj=jj+1;
        end
    end
    
    for session_number=session_numbers
        
        session_name=session_names{session_number};
        %% check if there is any analized video for for this experiment!
        dir_info=dir([BlDir,'\',sessions_strct(sessions_dirindx(session_number)).name,'\Record Node 101\neon_opto*.mat']);
        vid_exist=0;
        for vid_counter = 1: length(anlyized_video_list)
            if strcmp(anlyized_video_list(vid_counter).date(1:15),dir_info.date(1:15))
                disp(myKsDir)
                disp(anlyized_video_list(vid_counter).name)
                vid_exist = vid_exist + 1;
                vid_ID = vid_counter;
            end
        end
        if vid_exist
            %% load eye-tracking data and preprossesing
            liklihood_cut = .9;
%             vidObj = VideoReader(fullfile(Video_folder,anlyized_video_list(vid_ID).name));%'Y:\neon_pupil-Alireza-2022-09-01\videos\2020-04-14-18-21-34.mp4');
            T = readtable(fullfile(Video_folder,[anlyized_video_list(vid_ID).name(1:19),'DLC_effnet_b0_neon_pupilSep1shuffle1_500000_filtered.csv']),'PreserveVariableNames',true);%'Y:\neon_pupil-Alireza-2022-09-01\videos\2020-04-14-18-21-34DLC_effnet_b0_neon_pupilSep1shuffle1_500000.csv','PreserveVariableNames',true);
            T = table2array(T);
            %         T_copy = T;
%             ss = size (T);
            FPS = 45;
            %         FPS = vidObj.FrameRate;  %% video frame rate
%             if vidObj.NumFrames ~= ss(1)
%                 error(" missmatching frames in video and deeplabcut results")
%             else
%                 disp(' equal number of frame has been detected in both video file and and deeplabcut results')
%             end
            
            T(T(:,4) <liklihood_cut,2)  = nan; T(T(:,4) <liklihood_cut,3)  = nan;
            T(T(:,7) <liklihood_cut,5)  = nan; T(T(:,7) <liklihood_cut,6)  = nan;
            T(T(:,10)<liklihood_cut,8)  = nan; T(T(:,10)<liklihood_cut,9)  = nan;
            T(T(:,13)<liklihood_cut,11) = nan; T(T(:,13)<liklihood_cut,12) = nan;
            T(T(:,16)<liklihood_cut,14) = nan; T(T(:,16)<liklihood_cut,15) = nan;
            
            T(:,1) = T(:,1)./FPS; %% convert frame index to second
            
            h_diameter = sqrt((T(:,5)-T(:,8)).^2+(T(:,6)-T(:,9)).^2);
            v_diameter = sqrt((T(:,11)-T(:,14)).^2+(T(:,12)-T(:,15)).^2);
            
            h_diameter = smoothdata(h_diameter,'gaussian',5); % by default omits nans! Gaussian-weighted moving average over each window of A.
            v_diameter = smoothdata(v_diameter,'gaussian',5); % If a window contains all NaN values, then smoothdata returns NaN.
            diameters  = vertcat(h_diameter',v_diameter');
            ave_diameter = nanmean(diameters) ;
            
            x_center = (T(:,5)+T(:,8)+ T(:,11)+T(:,14))./4;
            y_center = (T(:,6)+T(:,9)+ T(:,12)+T(:,15))./4;
            center_loc = [x_center,y_center];
            center_loc = center_loc';
            x_nasal = nanmean(T(:,2));
            y_nasal = nanmean(T(:,3));
            mean_center_loc = nanmean(center_loc);
            Rot_ang = atand((mean_center_loc(2)-y_nasal)/(mean_center_loc(1)-x_nasal));
            rot_mat = [cosd(Rot_ang) -sind(Rot_ang); sind(Rot_ang) cosd(Rot_ang)]; %CREATE THE ROTATION MATRIX
            
            center_loc_rotatead = rot_mat*center_loc;
            
            %% stimate center and radius using circle fit
            center_loc_f = [];
            radius       = [];
            for fr_num = 1:length(T(:,1))
                XY = [T(fr_num,5:6);T(fr_num,8:9);T(fr_num,11:12);T(fr_num,14:15)];
                Par = CircleFitByPratt(XY);
                center_loc_f = horzcat(center_loc_f,Par(1:2)');
                radius       = horzcat(radius,Par(3));
            end
            ave_diameter_f = radius.*2;
            
            %% getting the video recording onsets and offset
            [d,timestamps, ~] = load_open_ephys_data_faster([BlDir,'\',sessions_strct(sessions_dirindx(session_number)).name,'\Record Node 101\100_ADC5.continuous']);
            fs = 30000;    % sampling frequency %Hz
            d_ind = (d  > max(d)/2);
            Video_trigger_ind = find(diff(d_ind));
            Video_trigger_time = timestamps(Video_trigger_ind); %#ok<FNDSB>
            T(:,1) = T(:,1)+ Video_trigger_time(1);
            
            %% geting the stimulus onsets
            if session_number==1 || session_number==3
                ITI       = .09;
                ver_thr   = 0.9; %vertical threshold (as a factor of std)
                hor_thr   = (ITI*fs);
            else
                ITI                  = 0.2;  %inter trial interval in second
                ver_thr              = 0.5; %vertical threshold (as a factor of std)
                hor_thr              = (ITI*fs);
            end
            [d,photo_timestamps, ~] = load_open_ephys_data_faster([BlDir,'\',sessions_strct(sessions_dirindx(session_number)).name,'\Record Node 101\100_ADC1.continuous']);
            d = double(d); d = smoothdata(d,'movmean',7);   % smooth data
            d = (d-mean(d))/std(d);     % Z_score data
            ind = find(d>ver_thr);      % apply threshold for value of data as a factor of std
            dif = diff(ind); dif(dif<hor_thr)=0;
            [~,locs] = findpeaks(dif);
            SOI = ind(locs+1); SOI = SOI';    % stimuli onset index
            SOI = [ind(1),SOI];               % #ok<AGROW>
            Stim_onset_time = photo_timestamps(SOI);
            
            %% trial schedule
            trial_schedule = stimulus_schedule{session_number};
            Nt  = min([length(trial_schedule),length(stimulus_time{session_number})-2]);   %#ok<*IDISVAR>
            pst = 0.2;                                                % pre-stimulus time [sec]
            stim_len = mean(stimulus_length{session_number})/fs;
            if strcmp(session_name,'neon')
                pst = 0.3;                                              % pre-stimulus time sec
                trial_schedule(Nt+1:end,:)=[];
                if min(size(trial_schedule))==2
                    realcolor_loc = trial_schedule(:,1)==2;
                    trial_schedule(realcolor_loc,2) = trial_schedule(realcolor_loc,2).*10;
                    trial_schedule = trial_schedule(:,2);
                elseif  min(size(trial_schedule))==3
                    phys_loc        = (trial_schedule(:,1)==1 & trial_schedule(:,2)==2);
                    ring_ctrl_loc   = (trial_schedule(:,1)==2 & trial_schedule(:,2)==1);
                    Square_ctrl_loc = (trial_schedule(:,1)==2 & trial_schedule(:,2)==2);
                    circle_ctrl_loc = (trial_schedule(:,1)==3 & trial_schedule(:,2)==1);
                    circle_phy_loc  = (trial_schedule(:,1)==3 & trial_schedule(:,2)==2);
                    
                    trial_schedule(phys_loc,3)        = trial_schedule(phys_loc,3).*10;
                    trial_schedule(ring_ctrl_loc,3)   = trial_schedule(ring_ctrl_loc,3).*100;
                    trial_schedule(Square_ctrl_loc,3) = trial_schedule(Square_ctrl_loc,3).*1000;
                    trial_schedule(circle_ctrl_loc,3) = trial_schedule(circle_ctrl_loc,3).*10000;
                    trial_schedule(circle_phy_loc,3)  = trial_schedule(circle_phy_loc,3).*100000;
                    
                    trial_schedule=trial_schedule(:,3);
                elseif  min(size(trial_schedule)) == 4
                    
                    phys_loc        = (trial_schedule(:,1)==1 & trial_schedule(:,2)==2 & trial_schedule(:,4)==0);
                    ring_ctrl_loc   = (trial_schedule(:,1)==2 & trial_schedule(:,2)==1);
                    Square_ctrl_loc = (trial_schedule(:,1)==2 & trial_schedule(:,2)==2);
                    circle_ctrl_loc = (trial_schedule(:,1)==3 & trial_schedule(:,2)==1);
                    circle_phy_loc  = (trial_schedule(:,1)==3 & trial_schedule(:,2)==2);
                    neon_opto_loc   = (trial_schedule(:,1)==1 & trial_schedule(:,2)==1 & trial_schedule(:,4)==1);
                    phys_opto_loc   = (trial_schedule(:,1)==1 & trial_schedule(:,2)==2 & trial_schedule(:,4)==1);
                    
                    trial_schedule(phys_loc,3)        = trial_schedule(phys_loc,3).*10;
                    trial_schedule(ring_ctrl_loc,3)   = trial_schedule(ring_ctrl_loc,3).*100;
                    trial_schedule(Square_ctrl_loc,3) = trial_schedule(Square_ctrl_loc,3).*1000;
                    trial_schedule(circle_ctrl_loc,3) = trial_schedule(circle_ctrl_loc,3).*10000;
                    trial_schedule(circle_phy_loc,3)  = trial_schedule(circle_phy_loc,3).*100000;
                    trial_schedule(neon_opto_loc,3)   = trial_schedule(neon_opto_loc,3).*1000000;
                    trial_schedule(phys_opto_loc,3)   = trial_schedule(phys_opto_loc,3).*10000000;
                    
                    trial_schedule = trial_schedule(:,3);
                end
            end
            %%
            
            
            condLabels = unique(trial_schedule);
            trial_start_time = Stim_onset_time-pst;        % starting time of trials in the data
            %         trial_end_time   = trial_start_time(2:end);    % ending time of trials.
            Stim_offset_time = Stim_onset_time+stim_len;   % stimulus offset
            trial_len = median(diff(trial_start_time));    % in second
            Time_ax = -pst:1/FPS:trial_len;
            trial_len_ind = length(Time_ax);
            Diameter_data    = [];
            Diameter_data_f  = [];
            pupil_diameter   = cell(2,1);
            pupil_diameter_f = [];
            pupil_center     = cell(2,1);
            pupil_center_f   = [];            
            pupil_center_rot = cell(2,1);
            tril_counter     = 0;
            %% loop over all condition to restructure the data
            for cond = 1:length(condLabels)
                ind_trGroup = trial_schedule(1:Nt)==condLabels(cond);
                cond_trial_onsets = trial_start_time(ind_trGroup);
                %cond_trial_offset = trial_end_time(ind_trGroup);
                cond_stim_onsets  = Stim_onset_time(ind_trGroup);
                cond_stim_offset  = Stim_offset_time(ind_trGroup);
                con_Ntrial(cond)  = length(cond_trial_onsets);
                for trial_number = 1:con_Ntrial(cond) %% loop over all trails in a condition for alignment
                    tril_counter  = tril_counter +1;
                    Trial_ind = T(:,1) >= cond_trial_onsets(trial_number) &  T(:,1) <= cond_trial_onsets(trial_number)+trial_len;
                    Trial_strat_ind = find( T(:,1) >= cond_trial_onsets(trial_number),1);
                    Stim_ind  = T(:,1) >= cond_stim_onsets(trial_number) &  T(:,1) <= cond_stim_offset(trial_number);
                    Blank_ind = T(:,1) >= cond_trial_onsets(trial_number) &  T(:,1) <= cond_stim_onsets(trial_number);
                    %% pupil info  calculation
                    Trial_pupil_size    = ave_diameter(Trial_strat_ind:Trial_strat_ind+trial_len_ind-1);
                    Stim_pupil_size     = ave_diameter(Stim_ind);
                    Blank_pupil_size    = ave_diameter(Blank_ind);
                    pupil_diameter{tril_counter} = Trial_pupil_size;
                    pupil_center{tril_counter} = center_loc(:,Trial_ind);
                    pupil_center_rot{tril_counter} = center_loc_rotatead(:,Trial_ind);
                    Diameter_data = [Diameter_data; cond, trial_number, nanmean(Trial_pupil_size), nanmean(Stim_pupil_size), nanmean(Blank_pupil_size),double(sum(isnan(Trial_pupil_size))>5)]; %#ok<AGROW>

                    Trial_pupil_size    = ave_diameter_f(Trial_strat_ind:Trial_strat_ind+trial_len_ind-1);
                    Stim_pupil_size     = ave_diameter_f(Stim_ind);
                    Blank_pupil_size    = ave_diameter_f(Blank_ind);
                    pupil_diameter_f    = [pupil_diameter_f;Trial_pupil_size];
                    pupil_center_f      = [pupil_center_f;center_loc_f(:,Trial_strat_ind:Trial_strat_ind+trial_len_ind-1)];
                    Diameter_data_f = [Diameter_data_f; cond, trial_number, nanmean(Trial_pupil_size), nanmean(Stim_pupil_size), nanmean(Blank_pupil_size),double(sum(isnan(Trial_pupil_size))>5)]; %#ok<AGROW>

                end
            end
            
            Pupil_data.diameter_data            = Diameter_data;
            Pupil_data.Diameter_data_f          = Diameter_data_f;
            Pupil_data.detail_of_diameter_data  = "colls : cond, trial number, ave-diameter (trial duration), ave-diameter (stim duration), ave-diameter (pre stim duration), nan flag";
            Pupil_data.pupil_diameter           = pupil_diameter;
            Pupil_data.pupil_diameter_f         = pupil_diameter_f;
            Pupil_data.detail_of_trial_diameter = "cells = different trials";
            Pupil_data.time_ax                  = Time_ax;
            Pupil_data.FrameRate                = FPS;
            Pupil_data.pupil_center             = pupil_center;
            Pupil_data.pupil_center_f             = pupil_center_f;
            Pupil_data.pupil_center_rot         = pupil_center_rot;

            save([myKsDir,'\',session_name,'Pupil_data.mat'],'Pupil_data')
            
            figure("Visible","off")
            histogram(Diameter_data(:,1)); xlabel('average disk speed')
            saveas(gcf,[myKsDir,'\',session_name,' Pupil_data_hist.jpg'])
            close gcf
        end
        
    end
end