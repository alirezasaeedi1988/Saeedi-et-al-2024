function f_RF_single_unit_neonFull(myKsDir,BlDir,Plot)
 
% myKsDir     = KS_dirs{day}%'D:\matwork\Data\higher_vision\M164_2020_02_12';
% BlDir       = BL_dirs{day}
% SD=2.5; %masking treshold
SD=3; %masking treshold

%% load Stimulus data and find different sessions!
stimTimeName = dir(fullfile(myKsDir, '*time.mat'));
if ~isempty(stimTimeName)
    stimTimeName = fullfile(myKsDir, stimTimeName(1).name);
    load(stimTimeName)
else
    error('there is no Stimulus-time in the directory');
end
schedule_length = cellfun(@length,stimulus_schedule,'uni',false); %#ok<*USENS>
schedule_length=cell2mat(schedule_length);
RF_sess=find(schedule_length>10000);
% size_sess=find(schedule_length==1200);

%%  laod channel map
chanMapName = dir(fullfile(myKsDir, 'chan*.mat'));
if ~isempty(chanMapName)
    chanMapName = fullfile(myKsDir, chanMapName(1).name); load(chanMapName)
else
    error('there is no channel map in the directory');
end



%% laod kilosort results
spike_times    = readNPY([myKsDir ,'\spike_times.npy']);
spike_cluster  = readNPY([myKsDir ,'\spike_clusters.npy']);
% [cids, cgs]    = readClusterGroupsCSV([myKsDir ,'\cluster_group.tsv']);

%% laod neon data
neonDataName = dir(fullfile(myKsDir, 'neon_spike_data.mat'));
if ~isempty(neonDataName)
    neonDataName = fullfile(myKsDir, neonDataName(1).name); load(neonDataName)
else
    error('there is no neon_data in the directory');
end
evoked_cids = cell2mat(spike_data.evoked_cids);

%% load log file for size modulation or random dot exp to get the monitor info
logFileName = dir(fullfile(BlDir,'Lsize*.mat'));
TF=dir(fullfile(BlDir,'random_dot*.mat'));

if ~isempty(logFileName)
    logData=load(fullfile(BlDir,logFileName(1).name));
    if isfield(logData,'distance_from_screen')
        monitor_dist=logData.distance_from_screen; %cm
    else
        monitor_dist=12;  %cm
    end
    monitor_X      = logData.screen_width_in_cm;%51.8; %cm
    monitor_Y      = 32.4; %cm
    
elseif ~isempty(TF)
    logFileName = dir(fullfile(BlDir,'random_dot*.mat'));
    logData=load(fullfile(BlDir,logFileName(1).name));
    if isfield(logData.scr,'dist')
        monitor_dist=logData.scr.dist; %cm
    else
        monitor_dist=12;  %cm
    end
    monitor_Y      = logData.scr.height;%32.4; %cm
    monitor_X      = logData.scr.width; %51.8; %cm
else  %% check open Ephys folder to find logdata
    dirinfo=dir(BlDir);
    dirinfo(~[dirinfo.isdir]) = [];
    for i=1:length(dirinfo)
        if length(dirinfo(i).name)>4
            mpath=dir(fullfile(BlDir,dirinfo(i).name,'Record Node 101','random_dot*.mat'));
            if ~isempty(mpath)
                matfile_path = fullfile(mpath.folder,mpath.name);
                logData=load(matfile_path);
                break
            end
        end
    end
    %     if isempty(mpath)
    %         error('no matfile in \n %s \n', BlDir)
    %     end
    if isfield(logData.scr,'dist')
        monitor_dist=logData.scr.dist; %cm
    else
        monitor_dist=12;  %cm
    end
    monitor_Y      = logData.scr.height;%32.4; %cm
    monitor_X      = logData.scr.width; %51.8; %cm
end

if ~exist('monitor_dist','var')
    monitor_dist = 12;  %cm
    monitor_X    = 51.8;
    monitor_Y      = 32.4; %cm
end



%% other params
res_scale      = 8;
min_lim_size   = 5;   % degree
% max_lim_size   = 25;  % degree


for sess=1:length(RF_sess)
% if length(RF_sess)==2
  
    session_name='RF';
    if sess==2
        session_name=[session_name,'2'];
    end
    outputDir           = fullfile(myKsDir,session_name);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir)
    end
    Stim_onset_ind = floor(stimulus_point{RF_sess(sess)});
    trial_schedule = stimulus_schedule{RF_sess(sess)};
    nsample        = floor((mean(diff(Stim_onset_ind(2:end-1))))/2);      %number of samples for each trials
    Nt             = min([length(trial_schedule),length(Stim_onset_ind)-2]);
    grid_cols      = max(unique(trial_schedule(:,1)));
    grid_rows      = max(unique(trial_schedule(:,2)));
    col_width      = monitor_X/(grid_cols*res_scale);  % cm
    row_width      = monitor_Y/(grid_rows*res_scale);  % cm
    col_angle      = 2*atand(col_width/(2*monitor_dist));
    row_angle      = 2*atand(row_width/(2*monitor_dist));
    
    
    
    %% trial_counter
    trial_counter  = zeros(grid_rows,grid_cols);
    for trial_idx  = 1:numel(Stim_onset_ind)
        trial_counter(trial_schedule(trial_idx,2), trial_schedule(trial_idx,1)) = ...
            trial_counter(trial_schedule(trial_idx,2), trial_schedule(trial_idx,1)) ...
            + 1;
    end
    
    
    %% ring polygons
    [~,ringimg_filled]    = f_ringimg();
    ringim_size=size(ringimg_filled);
    ringPixWidth  = monitor_X/max(ringim_size);
    ringPixHeight = monitor_Y/min(ringim_size);
    
    BW = imbinarize(ringimg_filled);
    BW_filled = imfill(BW,'holes');
    boundaries = bwboundaries(BW_filled);
    circles={};
    circles_cntr = [];
    for i=1:length(boundaries)
        circles{i}=polyshape(boundaries{i}(:,2)*ringPixWidth,boundaries{i}(:,1)*ringPixHeight); %#ok<*AGROW>
        [x,y] = centroid(circles{i});
        circles_cntr = [circles_cntr;[x,y]];
%         plot(circles{i}); hold on;  plot(x,y,'r*')

    end
    
    %% define empty matrices
    RF_size   = nan(length(evoked_cids),1);
    RSS       = nan(length(evoked_cids),1); % goodness of fitting (1-(Residual_sum_of_squares[fit]/Residual_sum_of_squares[tot]))
    RF_HWHM_X = nan(length(evoked_cids),1);
    RF_HWHM_Y = nan(length(evoked_cids),1);
    RF_loc    = nan(length(evoked_cids),2);
    RF_dist   = nan(length(evoked_cids),1);
    X_Bound   = {};%zeros(length(evoked_cids),1);
    Y_Bound   = {};%zeros(length(evoked_cids),1);
    
    RF_intrsct_ratio= nan(length(evoked_cids),1);
    RF_intrsct_ratio_ellips= nan(length(evoked_cids),1);
    
    RF_area   = nan(length(evoked_cids),1);
    RF_area_deg   = nan(length(evoked_cids),1);
    RF_area_ellips   = nan(length(evoked_cids),1);
    fs = 30000;
    dx=.001;
    stim_len = nsample/fs;
    xp = (-nsample/2)/fs:dx:(nsample*3/2)/fs;
    PSTH = zeros(length(evoked_cids),length(xp));
    spikes  = cell(1,length(evoked_cids));
    latency = nan(length(evoked_cids),1);
    
    [~, onset_ind]     = min(abs(xp));
    [~,offset_ind]     = min(abs(xp-stim_len));

    for neuron=1:length(evoked_cids) %% loop over clusters, spike time for selected cluster is extracted
        tr_count = zeros(grid_rows,grid_cols);
        spike = [];
        firnig_rate = zeros(grid_rows, grid_cols);
        cid         = spike_cluster==evoked_cids(neuron);
        Neuron_spike_times = spike_times(cid);
        Neuron_spike_times = Neuron_spike_times';
        Neuron_spike_times = double(Neuron_spike_times);
        for trial_number = 1:Nt
            tr_count(trial_schedule(trial_number,2), trial_schedule(trial_number,1)) = ...
                tr_count(trial_schedule(trial_number,2), trial_schedule(trial_number,1)) + 1;
            
            
            Trial_spike_ind   = Neuron_spike_times>= Stim_onset_ind(trial_number)-nsample/2 & ...
                Neuron_spike_times<Stim_onset_ind(trial_number)++(nsample*3/2);
            Trial_spike_times =  Neuron_spike_times(Trial_spike_ind);
            tr_N_spike = length(Trial_spike_times);
            if ~isempty(Trial_spike_times)
                %%% align
                Trial_spike_times=(Trial_spike_times-Stim_onset_ind(trial_number))./(fs);
                
                spike=[spike;...
                    Trial_spike_times' ,...
                    ones(tr_N_spike,1)*tr_count(trial_schedule(trial_number,2), trial_schedule(trial_number,1)) ,...
                    ones(tr_N_spike,1)*trial_schedule(trial_number,2),...
                    ones(tr_N_spike,1)*trial_schedule(trial_number,1)]; 
                
            end
            
            firnig_rate(trial_schedule(trial_number,2), trial_schedule(trial_number,1))=...
                firnig_rate(trial_schedule(trial_number,2), trial_schedule(trial_number,1))+ ...
                length(Trial_spike_times)/((trial_counter(trial_schedule(trial_number,2), trial_schedule(trial_number,1)))*(nsample/fs));
            
        end
        if ~std(firnig_rate(:))==0
            firnig_rate2 = (firnig_rate-mean(firnig_rate(:)))/std(firnig_rate(:));
            firnig_rate2(abs(firnig_rate2)<=2.25) = 0;
            [Rows,Cols] = find(firnig_rate2);
            RC = logical(firnig_rate2);
            hit_spikes=[]; %% spike in responses to hit squares!
            for r =1:length(Cols)
                idx = spike(:,3)== Rows(r) & spike(:,4)== Cols(r);
                hit_spikes = [hit_spikes; spike(idx,1)];
            end
            PSTH(neuron,:) = hist(hit_spikes,xp)./(sum(tr_count(RC))*dx);
            RelFreq = hist(hit_spikes,xp);
            RelFreq = movsum(RelFreq,2);
            spikes{neuron} = spike;
            teta= f_resp_latency(RelFreq,onset_ind,0,offset_ind);
            if ~isempty(teta)
                latency(neuron)=teta;
            end
            if Plot
                figure('Visible','off'); bar(xp,PSTH(neuron,:))
                xline(0,'-.g','LineWidth',1); xline(0.1,'-.g','LineWidth',1)
                if ~isempty(teta)
                    xline(latency(neuron)*dx,'-.r',[num2str(latency(neuron)),' ms'],'LineWidth',1,'LabelVerticalAlignment','top','LabelOrientation','horizontal');
                end
                ylabel('Firing Rate (Hz)');    ax = gca; ax.FontSize = 10;
                xlabel('Time (sec)')
                title(num2str(evoked_cids(neuron)))
                saveas(gcf,fullfile(outputDir,['RF_psth-',num2str(evoked_cids(neuron)),'.jpg']))
                saveas(gcf,fullfile(outputDir,['RF_psth-',num2str(evoked_cids(neuron)),'.pdf']))
                close gcf
            end
        end
        S = imresize(firnig_rate,res_scale);
        if ~std(S(:))==0
            S=(S-mean(S(:)))/std(S(:));S(abs(S)<= SD)=0; S=abs(S);
            if ~std(S(:))==0
                %% get the boundaries of RF
                BW = imbinarize(S);
                BW_filled = imfill(BW,'holes');
                boundaries2 = bwboundaries(BW_filled,4);
                
                
                
                %% 2D guassian fitting and finding the FR properties
                [A,RSS(neuron),Boundaries]=f_Gauss2D(S,Plot,0,outputDir,evoked_cids(neuron),row_angle ,col_angle,boundaries2,col_width,row_width);
                %[Amp,x0,wx,y0,wy,theta]=f_Gauss2D(S,Plot,FitOrientation,outputDir,cell_id,row_angle ,col_angle);
                %             if isempty(Boundaries); Boundaries=zeros(1,2);   end
                HWHM_X       = sqrt(2*log(2))*A(3);        % half width at half maximum
                HWHM_Y       = sqrt(2*log(2))*A(5);        % half width at half maximum
                HWHM_Y_angle = HWHM_Y*row_angle;
                HWHM_X_angle = HWHM_X*col_angle;
                RF_size(neuron)    = (HWHM_Y_angle+HWHM_X_angle)/2; %angle (degree)
                %RF_size(neuron)    = pi*HWHM_X*HWHM_Y;
                RF_loc(neuron,:)   = [A(2)*col_width, A(4)*row_width]; %cm
                RF_HWHM_X(neuron)  = HWHM_X*col_width; %cm
                RF_HWHM_Y(neuron)  = HWHM_Y*row_width; %cm
                RF_dist(neuron) = min(vecnorm(circles_cntr- RF_loc(neuron,:),1,2)); %% this the Euclidian distanc of RF center to closest inducer patch in cm.
                
                X_Bound{neuron}    = Boundaries(:,2)*col_width; %cm
                Y_Bound{neuron}    = Boundaries(:,1)*row_width; %cm
                %% area intersection
                ppp=polyshape(X_Bound{neuron},Y_Bound{neuron} );
                RF_area(neuron) = area(ppp);  % cm^2
                area_int=0;
                for cir=1:length(circles)
                    polyout = intersect(ppp,circles{cir});
                    area_int= area_int + area(polyout);
                end
                RF_intrsct_ratio(neuron) = area_int/RF_area(neuron)*100;
                %% RF area size in degree
                X_Bound_ang    = Boundaries(:,2)*col_angle; %degree
                Y_Bound_ang    = Boundaries(:,1)*row_angle; %degree
                ppp_ang=polyshape(X_Bound_ang,Y_Bound_ang);
                RF_area_deg(neuron) = area(ppp_ang);  % degree^2
                %% ellips
                t = linspace(0,2*pi) ;
                ellips_x = A(2)+ HWHM_X*cos(t) ;
                ellips_y = A(4)+ HWHM_Y*sin(t) ;
                ellips_x = ellips_x*col_width; %cm
                ellips_y = ellips_y*row_width;%cm
                %% area intersection for ellips
                ppp=polyshape(ellips_x,ellips_y );
                RF_area_ellips(neuron) = area(ppp);  % cm^2
                area_int=0;
                for cir=1:length(circles)
                    polyout = intersect(ppp,circles{cir});
                    area_int= area_int + area(polyout);
                end
                RF_intrsct_ratio_ellips(neuron) = area_int/RF_area_ellips(neuron)*100;
                
            end
        end
    end
    RF_size(RF_size<min_lim_size) = nan;
%     RF_size(RF_size>max_lim_size) = nan;
    % good_rf = ~isnan(RF_size);
    RF_loc_shifted      = RF_loc;
    RF_loc_shifted(:,1) = RF_loc_shifted(:,1)- grid_cols;
    RF_loc_shifted(:,2) = RF_loc_shifted(:,2)- grid_rows;
    
    
    RF_data.distance   = RF_dist;
    RF_data.spikes        = spikes;
    RF_data.PSTH          = PSTH;
    RF_data.xp            = xp;
    RF_data.dx            = dx;
    RF_data.stim_len      = stim_len;
    RF_data.latency       = latency;
    RF_data.intrsct_ratio = RF_intrsct_ratio;
    RF_data.intrsct_ratio_ellips = RF_intrsct_ratio_ellips;
    RF_data.area         = RF_area;
    RF_data.area_deg     = RF_area_deg;
    RF_data.area_ellips  = RF_area_ellips;
    RF_data.HWHM_X       = RF_HWHM_X;
    RF_data.HWHM_Y       = RF_HWHM_Y;
    RF_data.X_Bound      = X_Bound;
    RF_data.Y_Bound      = Y_Bound;
    RF_data.RF_size      = RF_size;
    RF_data.RF_loc       = RF_loc;
    RF_data.RF_loc_shift = RF_loc_shifted;
    RF_data.RSS          = RSS;
    RF_data.res_scale    = res_scale;
    RF_data.session_name = session_name;
    RF_data.evoked_cids  = evoked_cids;
    RF_data.monitor_size = [monitor_X, monitor_Y]; %%cm
    save(fullfile(myKsDir,[session_name,'_data.mat']),'RF_data')
end
