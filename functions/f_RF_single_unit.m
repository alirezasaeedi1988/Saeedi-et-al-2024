function f_RF_single_unit(myKsDir,BlDir,session_number,session_name,Plot)


% myKsDir     = 'D:\matwork\Data\higher_vision\M164_2020_02_12';
% session_number = 1;
% session_name='RF';

outputDir           = fullfile(myKsDir,session_name);
    if ~exist(outputDir, 'dir')
       mkdir(outputDir)
    end

%%  laod channel map
chanMapName = dir(fullfile(myKsDir, 'chan*.mat'));
if ~isempty(chanMapName)
    chanMapName = fullfile(myKsDir, chanMapName(1).name); load(chanMapName)
else
    error('there is no channel map in the directory');
end

%% load Stimulus data
stimTimeName = dir(fullfile(myKsDir, '*time.mat'));
if ~isempty(stimTimeName)
    stimTimeName = fullfile(myKsDir, stimTimeName(1).name); load(stimTimeName)
else
    error('there is no Stimulus-time in the directory');
end

%% laod kilosort results
spike_times    = readNPY([myKsDir ,'\spike_times.npy']);
spike_cluster  = readNPY([myKsDir ,'\spike_clusters.npy']);
% [cids, cgs]    = readClusterGroupsCSV([myKsDir ,'\cluster_group.tsv']);

%% laod neon data
neonDataName = dir(fullfile(myKsDir, 'neon*.mat'));
if ~isempty(neonDataName)
    neonDataName = fullfile(myKsDir, neonDataName(1).name); load(neonDataName)
else
    error('there is no neon_data in the directory');
end
evoked_cids = cell2mat(spike_data.evoked_cids);

%% load size modulation data to get the distance
logFileName = dir(fullfile(BlDir,'Lsize*.mat'));
if ~isempty(logFileName)
    logData=load(fullfile(BlDir,logFileName(1).name));
    if isfield(logData,'distance_from_screen')
        monitor_dist=logData.distance_from_screen; %cm
    end
    monitor_X      = logData.screen_width_in_cm;%51.8; %cm

else
    monitor_dist = 12;  %cm
    monitor_X    = 51.8;
end
%% other params
res_scale      = 8;
monitor_Y      = 32.4; %cm

min_lim_size   = 5;   % degree      
max_lim_size   = 25;  % degree

Stim_onset_ind = floor(stimulus_point{session_number});      
trial_schedule = stimulus_schedule{session_number};
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
for i=1:length(boundaries)
    circles{i}=polyshape(boundaries{i}(:,2)*ringPixWidth,boundaries{i}(:,1)*ringPixHeight); %#ok<*AGROW>
end

%% define empty matrices
RF_size   = zeros(length(evoked_cids),1);
RSS       = zeros(length(evoked_cids),1); % goodness of fitting (1-(Residual_sum_of_squares[fit]/Residual_sum_of_squares[tot]))
RF_HWHM_X = zeros(length(evoked_cids),1);
RF_HWHM_Y = zeros(length(evoked_cids),1);
RF_loc    = zeros(length(evoked_cids),2);
X_Bound   = {};%zeros(length(evoked_cids),1);
Y_Bound   = {};%zeros(length(evoked_cids),1);
RF_intrsct_ratio= zeros(length(evoked_cids),1);
RF_area   = zeros(length(evoked_cids),1);

 for neuron=1:length(evoked_cids) %% loop over clusters, spike time for selected cluster is extracted

     firnig_rate = zeros(grid_rows, grid_cols);
     cid         = spike_cluster==evoked_cids(neuron);
     Neuron_spike_times = spike_times(cid);
     Neuron_spike_times = Neuron_spike_times';
     Neuron_spike_times = double(Neuron_spike_times);
     for trial_number = 1:Nt
         Trial_spike_ind   = Stim_onset_ind(trial_number)<=Neuron_spike_times & ...
              Neuron_spike_times<Stim_onset_ind(trial_number)+nsample;
         Trial_spike_times =  Neuron_spike_times(Trial_spike_ind);                    
         firnig_rate(trial_schedule(trial_number,2), trial_schedule(trial_number,1))=...
             firnig_rate(trial_schedule(trial_number,2), trial_schedule(trial_number,1))+ ...
             length(Trial_spike_times)/((trial_counter(trial_schedule(trial_number,2), trial_schedule(trial_number,1)))*(nsample/fs));

     end
     
     S = imresize(firnig_rate,res_scale);
     if ~std(S(:))==0
         S=(S-mean(S(:)))/std(S(:));S(abs(S)<= 1.8)=0; S=abs(S);
		 
         %% get the boundaries of RF
         BW = imbinarize(S);
         BW_filled = imfill(BW,'holes');
         boundaries2 = bwboundaries(BW_filled,4);
         allLengths = cellfun(@length,boundaries2);
         [~,maxLoc] = max(allLengths);
         Boundaries = boundaries2{maxLoc};
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
         
         %% 2D guassian fitting and finding the FR properties
         [A,RSS(neuron)]=f_Gauss2D(S,Plot,0,outputDir,evoked_cids(neuron),row_angle ,col_angle,Boundaries );
         %[Amp,x0,wx,y0,wy,theta]=f_Gauss2D(S,Plot,FitOrientation,outputDir,cell_id,row_angle ,col_angle);
         HWHM_X       = sqrt(2*log(2))*A(3);        % half width at half maximum
         HWHM_Y       = sqrt(2*log(2))*A(5);        % half width at half maximum
         HWHM_Y_angle = HWHM_Y*row_angle;
         HWHM_X_angle = HWHM_X*col_angle;
         RF_size(neuron)    = (HWHM_Y_angle+HWHM_X_angle)/2; %angle (degree)
         %RF_size(neuron)    = pi*HWHM_X*HWHM_Y;
         RF_loc(neuron,:)   = [A(2)*col_width, A(4)*row_width]; %cm
         RF_HWHM_X(neuron)  = HWHM_X*col_width; %cm
         RF_HWHM_Y(neuron)  = HWHM_Y*row_width; %cm

     end
 end
RF_size(RF_size<min_lim_size) = nan;
RF_size(RF_size>max_lim_size) = nan;
% good_rf = ~isnan(RF_size);
RF_loc_shifted      = RF_loc;
RF_loc_shifted(:,1) = RF_loc_shifted(:,1)- grid_cols;
RF_loc_shifted(:,2) = RF_loc_shifted(:,2)- grid_rows;

RF_data.intrsct_ratio= RF_intrsct_ratio;
RF_data.area         = RF_area;
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

