function [preferred_angle,preferred_ind]=f_find_preferred6(Pdf,stim_angle,angle_loc,onset_ind,offset_ind,type)
%% find_preferred, inputs:Pdf,stim_angle,angle_loc,onset_ind,offset_ind
% type= 'all' , 'neon' , 'phy'
ave_activity_angle    = zeros(1,length(stim_angle));
pre=mean(mean(Pdf(:,1:onset_ind-1)));

for ang=1:length(stim_angle) %%% loop over angles to find the preferred angle
    switch type
        case 'all'
            p=mean(Pdf(angle_loc(ang,1:end),:)); %-1 is to reject ctrl
        case 'neon'
            p=(Pdf(angle_loc(ang,1),:)); 
        case 'phy'
            p=mean(Pdf(angle_loc(ang,[2,6]),:)); %-1 is to reject ctrl 
        case 'grat'
            p=(Pdf(angle_loc(ang,1),:)); 
        case 'rand'
            p=(Pdf(angle_loc(ang,2),:));             
    end

    %             pre=mean(p(1:onset_ind));
    p=abs((p-pre));
    ave_activity_angle(ang)=mean(p(onset_ind:offset_ind));
end
[~,preferred_ind]=max(ave_activity_angle);
preferred_angle = stim_angle(preferred_ind);