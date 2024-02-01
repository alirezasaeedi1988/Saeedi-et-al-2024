function [preferred_angle,preferred_ind,osi,dsi,Rsq1,osi2,dsi2,Rsq2,osi3,dsi3,ave_FR_angles]=...
    f_find_preferred6_cond(spike_data,stim_angle,angle_loc,onset_ind,offset_ind,cond,myKsDir,neuron,Report)
%% find_preferred,OSI and DSI
%% inputs  : spike_data, stim_angle, angle_loc,onset_ind,offset_ind,cond,myKsDir,neuron,cid,Report
%% outputs : preferred_angle,preferred_ind,osi,dsi
% cond= 0 for overall condition



id=spike_data.evoked_cids{neuron};
Pdf=spike_data.pdfs{neuron};

ave_activity_angle    = zeros(1,length(stim_angle));
ave_activity_angle2    = zeros(1,length(stim_angle));

pre=mean(mean(Pdf(:,1:onset_ind-1)));

for ang=1:length(stim_angle) %%% loop over angles to find the preferred angle
    switch cond 
        case 0 %zero for all
            p = mean(Pdf(angle_loc(ang,:),onset_ind:offset_ind));
        otherwise
            p = Pdf(angle_loc(ang,cond),onset_ind:offset_ind);
    end
    ave_activity_angle2(ang)= mean(p);
    p = p-pre;
    ave_activity_angle(ang)= mean(p);
end
ave_FR_angles = mean(ave_activity_angle);
[~,preferred_ind] = max(ave_activity_angle);
preferred_angle   = stim_angle(preferred_ind);

inv_activity_angle = ave_activity_angle.*sign(mean(ave_activity_angle)); %% to invert the activity for inhibited cells

%% OSI and DSI
x_fit   = 1:length(stim_angle);
y_fit   = circshift(inv_activity_angle,5-preferred_ind); %% shift the peak to the 5th element 
g       = @(A,X) A(1)*exp( -((X-5).^2/(2*A(2)^2)))+ A(3)*exp( -((X-5+4).^2/(2*A(4)^2)))+A(5); %% double gaussian with pre-defined mean valeus
A0      = [.5,1,.5,1, mean(y_fit)];  %% initial guess
[A,RSS] = lsqcurvefit(g,A0,x_fit,y_fit); %% fitting
Rsq1     = 1-(RSS/(rssq(inv_activity_angle)^2));
new_y   = A(1)*exp( -((x_fit-5).^2/(2*A(2)^2)))+ A(3)*exp( -((x_fit-5+4).^2/(2*A(4)^2)))+A(5); %% generat data with the fitted function
new_y   = circshift(new_y,-(5-preferred_ind)); %% shift the peak to original location

% ffit = fit(stim_angle',ave_activity_angle','gauss2');


opp_ind = mod(preferred_ind+4,8);
if opp_ind ==0
    opp_ind=8;
end
dsi     = (new_y(preferred_ind)-new_y(opp_ind))/(new_y(preferred_ind)+new_y(opp_ind));
% ort_ind = mod(preferred_ind+2,8);
% if ort_ind ==0
%     ort_ind=8;
% end
% osi     = (new_y(preferred_ind)-new_y(ort_ind))/(new_y(preferred_ind)+new_y(ort_ind));
shifted_y=circshift(new_y,1-preferred_ind); %% shift the peak to the 1th element 
shifted_y=mean(reshape(shifted_y,4,2),2); % reshape the activities to get the mean value for similar orientations!
osi     = (new_y(preferred_ind)-shifted_y(2))/(new_y(preferred_ind)+shifted_y(2));

%% OSI2 and DSI2
[~,preferred_ind2] = max(ave_activity_angle2);
x_fit2   = 1:length(stim_angle);
y_fit2   = circshift(ave_activity_angle2,5-preferred_ind2); %% shift the peak to the 5th element 
g       = @(A,X) A(1)*exp( -((X-5).^2/(2*A(2)^2)))+ A(3)*exp( -((X-5+4).^2/(2*A(4)^2)))+A(5); %% double gaussian with pre-defined mean valeus
A0      = [.5,1,.5,1, mean(y_fit2)];  %% initial guess
[A,RSS] = lsqcurvefit(g,A0,x_fit2,y_fit2); %% fitting
Rsq2     = 1-(RSS/(rssq(ave_activity_angle2)^2));
new_y2   = A(1)*exp( -((x_fit2-5).^2/(2*A(2)^2)))+ A(3)*exp( -((x_fit2-5+4).^2/(2*A(4)^2)))+A(5); %% generat data with the fitted function
new_y2   = circshift(new_y2,-(5-preferred_ind2)); %% shift the peak to original location

opp_ind = mod(preferred_ind2+4,8);
if opp_ind ==0
    opp_ind=8;
end
dsi2     = (new_y2(preferred_ind2)-new_y2(opp_ind))/(new_y2(preferred_ind2)+new_y2(opp_ind));
% ort_ind = mod(preferred_ind2+2,8);
% if ort_ind ==0
%     ort_ind=8;
% end
% osi2     = (new_y2(preferred_ind2)-new_y2(ort_ind))/(new_y2(preferred_ind2)+new_y2(ort_ind));

shifted_y=circshift(new_y2,1-preferred_ind2); %% shift the peak to the 1th element 
shifted_y=mean(reshape(shifted_y,4,2),2); % reshape the activities to get the mean value for similar orientations!
osi2     = (new_y2(preferred_ind2)-shifted_y(3))/(new_y2(preferred_ind2)+shifted_y(3));

%% OSI3 and DSI3
dsi3=sqrt(sum(ave_activity_angle2.*sind(stim_angle))^2+sum(ave_activity_angle2.*cosd(stim_angle))^2)/sum(ave_activity_angle2);
osi3=sqrt(sum(ave_activity_angle2.*sind(2.*stim_angle))^2+sum(ave_activity_angle2.*cosd(2.*stim_angle))^2)/sum(ave_activity_angle2);

%% plotting
if_cond= cond==1 | cond==2;
if Report && if_cond 
    session_name='neon';
    outputDir           = fullfile(myKsDir,session_name);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir)
    end
    figure('units','normalized','outerposition',[0 0 0.4 1],'Visible','off');
    subplot(3,1,1)
    plot(stim_angle,ave_activity_angle,'-k*','LineWidth', 2)
    hold on
    plot(stim_angle,inv_activity_angle,'-b*','LineWidth', 2)
    plot(stim_angle,new_y,'r','LineWidth', 2)
    legend('reduced FR','inverted FR','gauss2 fit','Location','northeastoutside')
    xticks(stim_angle)
    title(['neuron-',num2str(id),' cond=',num2str(cond),' RSQ=',num2str(Rsq1)])
    ax = gca;ax.FontSize = 14; ax.FontName="Arial";
    
    subplot(3,1,2)
    plot(stim_angle,ave_activity_angle2,'-k*','LineWidth', 2)
    hold on
    plot(stim_angle,new_y2,'r','LineWidth', 2)
    legend('raw FR','gauss2 fit','Location','northeastoutside')
    xticks(stim_angle)
    title(['neuron-',num2str(id),' cond=',num2str(cond),' RSQ=',num2str(Rsq2)])
    ax = gca;
    ax.FontSize = 12;
    ax.FontName="Arial";    
    subplot(3,1,3)
    polarplot(deg2rad([stim_angle,0]),[ave_activity_angle2,ave_activity_angle2(1)],'-k','LineWidth', 2)
    hold on
    polarplot(deg2rad([stim_angle,0]),[inv_activity_angle,inv_activity_angle(1)],'-b','LineWidth', 2)
    thetaticks(stim_angle)
    legend('raw FR','reduced inverted FR','Location','northeastoutside')

    title({['DSI1=',num2str(dsi,'%.2f'),'  OSI1=',num2str(osi,'%.2f')];...
        ['DSI2=',num2str(dsi2,'%.2f'),'  OSI2=',num2str(osi2,'%.2f')];...
        ['DSI3=',num2str(dsi3,'%.2f'),'  OSI3=',num2str(osi3,'%.2f')]})
    pax = gca;
    pax.FontSize = 12;
    pax.FontName="Arial";
    saveas(gcf,[outputDir,'\tunning_fit_',num2str(id),'_',num2str(cond),'.png'])
    close gcf
end
%% interpolation for activities
% interpolated_angle = min(stim_angle):5:max(stim_angle);
%
% interp_act = interp1(stim_angle,ave_activity_angle,interpolated_angle,'spline');
% [~,inter_pref_ind]=max(interp_act);
% inter_pref_angle = interpolated_angle(inter_pref_ind);
%
% g = @(A,X) A(1)*exp( -((X-inter_pref_angle).^2/(2*A(2)^2)))+ A(3)*exp( -((X-inter_pref_angle+180).^2/(2*A(4)^2)))+A(5);
% A0=[1,1,1,1, mean(ave_activity_angle)];
% [A,RSS] = lsqnonlin(g,A0,interpolated_angle,interp_act);
% Rsq=1-(RSS/(rssq(interp_act)^2));
% new_y=A(1)*exp( -((interpolated_angle-inter_pref_angle).^2/(2*A(2)^2)))+ A(3)*exp( -((interpolated_angle-inter_pref_angle+180).^2/(2*A(4)^2)))+A(5);
%
% plot(interpolated_angle,interp_act);xticks(stim_angle)
% hold on
% plot(interpolated_angle,new_y);xticks(stim_angle)