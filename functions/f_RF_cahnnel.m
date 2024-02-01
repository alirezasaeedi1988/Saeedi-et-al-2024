function f_RF_cahnnel(myKsDir,base_folder)
addpath('W:\projects\alireza\core_toolbox')
% base_folder = 'Y:\MPI_EphysData_2020\higher_vision\M164_left\M164_2020_02_13';
% myKsDir='D:\matwork\data\higher_vision\M164_2020_02_13';

TF=dir(fullfile(base_folder,'*.ns3'));
if ~isempty(TF)
    fprintf('ns3 data has been found in \n %s \n',base_folder)
    matfile_path = dir(fullfile(base_folder, 'retinotopy*.mat'));
    matfile_path = fullfile(matfile_path.folder,matfile_path.name);
    pd_path = fullfile(base_folder, 'datafile001.ns3');
    nd_path = fullfile(base_folder, 'datafile001.ns6');
else
    dirinfo=dir(base_folder);
    dirinfo(~[dirinfo.isdir]) = [];
    for i=1:length(dirinfo)
        if length(dirinfo(i).name)>4
            mpath=dir(fullfile(base_folder,dirinfo(i).name,'Record Node 101','retinotopy*.mat'));
            if isempty(mpath)
                error('no matfile in \n %s \n', base_folder)
            else
                matfile_path = fullfile(mpath.folder,mpath.name);
                pd_path = fullfile(mpath.folder, '100_ADC1.continuous');
                nd_path = fullfile(mpath.folder, '100_CH1.continuous');
                break
            end
            
        end
    end
end
% is there a channel map file in this folder?
chanmapfile = dir(fullfile(myKsDir, 'chan*.mat'));
if ~isempty(chanmapfile)
    chanmapfile = fullfile(myKsDir, chanmapfile(1).name);
    load(chanmapfile)
else
    error('there is no channelMap in MyKsDir')
end

min_ISI_seconds =.05;%0.1;
event_offset_seconds = 0.02;
pre_onset_seconds = 0.0;
post_onset_seconds = 0.15;
pd_channel = 65; % only applies to ns3
res_scale      = 2;
monitor_Y      = 32.4; %cm
monitor_X      = 51.8; %cm
monitor_dist   = 18;   %cm

%% get onsets
[event_onset_seconds,trial_miss] = get_analog_signal_onsets_170221(pd_path,pd_channel, min_ISI_seconds,0,0);
% if trial_miss
%     event_onset_seconds(trial_miss)=[];
% end
event_onset_seconds = event_onset_seconds + event_offset_seconds;





load(matfile_path)
grid_cols      = max(unique(trial_schedule(:,1)));
grid_rows      = max(unique(trial_schedule(:,2)));
col_width      = monitor_X/(grid_cols*res_scale);  % cm
row_width      = monitor_Y/(grid_rows*res_scale);  % cm
col_angle      = 2*atand(col_width/(2*monitor_dist));
row_angle      = 2*atand(row_width/(2*monitor_dist));


trial_counter = zeros(grid_cols,grid_rows);
for trial_idx = 1:numel(event_onset_seconds)
    trial_counter(trial_schedule(trial_idx,1), trial_schedule(trial_idx,2)) = ...
        trial_counter(trial_schedule(trial_idx,1), trial_schedule(trial_idx,2)) ...
        + 1;
end

max_num_trials = max(trial_counter(:));


%% filter params
order = 2;
fcutlow_band  = 1000;
fcutlow_low  = 140;
fcuthigh_band = 10000;
fs = 30000;
[b_band,a_band] = butter(order,[fcutlow_band,fcuthigh_band]/(fs/2), 'bandpass');
[b_low,a_low] = butter(order,fcutlow_low/(fs/2), 'low');


%% other params
ephys_sampling_freq = 30000; % can also be retrieved from file but whatever
pre_onset_samples = pre_onset_seconds * ephys_sampling_freq;
post_onset_samples = post_onset_seconds * ephys_sampling_freq;
num_samples = pre_onset_samples + post_onset_samples + 1;

%% main
if ~isempty(TF)
    ns6_file = openNSx('read',nd_path,'e:1:64');
end
Nchan=64;
for channel_idx = 1:Nchan
    if ~isempty(TF)
        data = ns6_file.Data(channel_idx,:);
    else
        data = load_open_ephys_data_faster([mpath.folder '\' sprintf('100_CH%d.continuous', channel_idx)]);
    end
    continuous_MUA = filtfilt(b_low, a_low, ...
        abs(filtfilt(b_band,a_band,data)));
    
    MUA = nan(grid_cols, grid_rows, max_num_trials, num_samples);
    
    trial_counter = zeros(grid_cols,grid_rows);
    
    for trial_idx = 1:numel(event_onset_seconds)-1
        trial_counter(trial_schedule(trial_idx,1), trial_schedule(trial_idx,2)) = ...
            trial_counter(trial_schedule(trial_idx,1), trial_schedule(trial_idx,2)) ...
            + 1;
        
        condition_trial_idx =  trial_counter(trial_schedule(trial_idx,1), trial_schedule(trial_idx,2));
        
        event_onset_sample = event_onset_seconds(trial_idx) * ephys_sampling_freq;
        
        continuous_chunk = continuous_MUA( ...
            round(event_onset_sample - pre_onset_samples): ...
            round(event_onset_sample + post_onset_samples));
        
        MUA(trial_schedule(trial_idx,1), trial_schedule(trial_idx,2), condition_trial_idx, :) = ...
            continuous_chunk;
        % size(MUA)
    end
    
%     min_MUA = inf;
%     max_MUA = -inf;
%     
    %figure
    for i = 1:grid_cols
        for j = 1:grid_rows
            condition_all_trials = squeeze(MUA(i,j,:,:));
            mean_MUA(i,j) = mean(nanmedian(condition_all_trials));
        end
    end
    
    all_mean_MUA(channel_idx,:,:) = mean_MUA;
    
end

RF_size   = zeros(Nchan,1);
RSS       = zeros(Nchan,1); % Residual sum of squares

RF_loc    = zeros(Nchan,2);

for channel_idx = 1:Nchan
    
    ch=chanMap(channel_idx); %%top to down
    MUA_cha = squeeze(all_mean_MUA(ch,:,:))';
    
    S = imresize(MUA_cha,res_scale);
    if ~std(S(:))==0
        S=(S-mean(S(:)))/std(S(:));S(abs(S)<= 1.8)=0; S=abs(S);
        %smoothed_fr=imgaussfilt(firnig_rate);      %smoothing with gaussian filter
        [A,RSS(channel_idx)]=f_Gauss2D(S,0,0,base_folder,ch,row_angle ,col_angle );
        %[Amp,x0,wx,y0,wy,theta]=f_Gauss2D(S,Plot,FitOrientation,outputDir,cell_id,row_angle ,col_angle);
        HWHM_X       = sqrt(2*log(2))*A(3);        % half width at half maximum
        HWHM_Y       = sqrt(2*log(2))*A(5);        % half width at half maximum
        HWHM_Y_angle = HWHM_Y*row_angle;
        HWHM_X_angle = HWHM_X*col_angle;
        RF_size(channel_idx)    = (HWHM_Y_angle+HWHM_X_angle)/2;
        %RF_size(neuron)    = pi*HWHM_X*HWHM_Y;
        RF_loc(channel_idx,:)   = [A(2), A(4)];
    end
end

RF_loc_shifted      = RF_loc;
RF_loc_shifted(:,1) = RF_loc_shifted(:,1)- grid_cols;
RF_loc_shifted(:,2) = RF_loc_shifted(:,2)- grid_rows;


%% saving
RF_chan.RF_size      = RF_size;
RF_chan.RF_loc       = RF_loc;
RF_chan.RF_loc_shift = RF_loc_shifted;
RF_chan.RSS          = RSS;
RF_chan.res_scale    = res_scale;
RF_chan.chanDepth    = ycoords;


%% ploting
chan_label = 1:length(RF_loc);chan_label=chan_label';
bad_fit_ind = RSS > 100;
disp(myKsDir)
hem=input('which hem.? (r/l) \n');
RF_loc(bad_fit_ind,:)   = [];
ycoords(bad_fit_ind)    = [];
chan_label(bad_fit_ind) = [];
y_loc = grid_rows*res_scale+1-RF_loc(:,2);
figure('units','normalized','outerposition',[0 0 1 1])
plot(RF_loc(:,1),y_loc,'-b*')
hold on
rectangle('Position',[1 1 25 15],'EdgeColor','r')
if strcmp(hem,'l')
    xlabel('nasal   <<<<------------------------->>>>   temporal')
elseif strcmp(hem,'r')
    xlabel('temporal   <<<<------------------------->>>>   nasal')
else
    disp('not defined hemisphere')
end
b = num2str(chan_label); c = cellstr(b);
dx = 0.1; dy = 0.3;
text(RF_loc(:,1)+dx, y_loc+dy, c);
xticks([])
yticks([])
xlim ([0,27])
ylim ([0,17])

saveas(gcf,fullfile(myKsDir,'Channel RF location.png'))



BB = smooth(RF_loc(:,1));%,'rlowess');
figure('units','normalized','outerposition',[0 0 1 1])
plot(RF_loc(:,1),ycoords,'-b*')
hold on;plot(BB,ycoords,'-r')
ylabel('channel depth (\mum)')
xlabel('horizontal location on screen')
text(BB+dx*2, ycoords+dy, c);
if strcmp(hem,'l')
    title('nasal   <<<<------------------------->>>>   temporal')
elseif strcmp(hem,'r')
    title('temporal   <<<<------------------------->>>>   nasal')
end
legend('data','smoothed data')

saveas(gcf,fullfile(myKsDir,'Channel RF_Horizontal location.png'))

%% reversal points
x=abs(RF_loc_shifted(:,1));
X = smooth(x,'rlowess');
% figure
% plot(X)
% hold on
% plot(x)
% plot(RF_loc_shifted(:,1))
[~,reversals] = findpeaks(X,'MinPeakDistance',14);
% findpeaks(X,'MinPeakDistance',14)
% text(reversals+.02,pks+2,num2str(reversals))
% ylabel('xlocation')
% xlabel('channel')
fprintf('reversals points are \n')
disp(reversals)
inp= input('enter new reversals or 0 to keep the old one \n');
if inp~=0
    reversals=inp;
end
RF_chan.reversals    = reversals;
save(fullfile(myKsDir,'channel_RF_data.mat'),'RF_chan')
close all
