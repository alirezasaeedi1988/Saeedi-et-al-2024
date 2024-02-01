function teta=f_resp_latency(y,onset_ind,transient_nbins,offset_ind)


% y                =  sum(RFs{8});
% [~, onset_ind]   =  min(abs(xp));
% [~,offset_ind]   =  min(abs(xp-stim_len));
lambda           =  mean(y(1:onset_ind));
% pre              =  sum (y(1:onset_ind));
% post             =  sum (y(onset_ind+1:onset_ind*2));
% if post<pre
% %     y       = abs(y-pre);
%     lambda  =  lambda/2;
% end
prob = poisspdf(y(onset_ind+transient_nbins:offset_ind),lambda);
sadom_loc=find(prob<.01);
cons_sadom_loc=sadom_loc(diff(sadom_loc)==1);
teta=cons_sadom_loc(find(prob(cons_sadom_loc(1:end-2)+2)<=0.05,1));