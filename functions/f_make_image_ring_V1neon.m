function [fullneon,filled_fullneon]= f_make_image_ring_V1neon(BlDir,neon_sess)

% output: fullneon is an image with the exact size (oixelwise) of monitor and the presented stimulus on it.
%         screen_size_in_pixels 




% make circular rings
%% load neon log file
mat_fnames=dir(fullfile(BlDir,'neon_color_spread*.mat'));
if ~isempty(mat_fnames)
    neon_exp = load(fullfile(mat_fnames.folder,mat_fnames.name));
    
else
    root_content = struct2table(dir(BlDir));
    root_content = root_content(root_content.isdir==true & cellfun(@length,root_content.name)>10,:);

    neon_path = fullfile(char(root_content.folder(neon_sess)),char(root_content.name(neon_sess)));
    mat_fnames=dir(fullfile(neon_path,'neon_color_spread*.mat'));

    neon_exp  = load(fullfile(mat_fnames.folder,mat_fnames.name));

end

ringwidth=neon_exp.ringwidth;  % if possible run 0.3 too
% wrect= neon_exp.wrect;
% screen_size_in_pixels = wrect(3:4);
screen_width_in_pixels  = neon_exp.wrect(3);
screen_height_in_pixels = neon_exp.wrect(4);
stimrect=neon_exp.stimrect;
horiz_shift=neon_exp.horiz_shift;
vert_shift=neon_exp.vert_shift;
new_stimrect = stimrect+floor([horiz_shift vert_shift horiz_shift vert_shift]);




ringsize=neon_exp.ringsize; 
rings=neon_exp.rings(3);

[x,y]=meshgrid(-ringsize:ringsize,-ringsize:ringsize);
[~,r]=cart2pol(x,y);
ringimg=x*0;

ringimg(rings-ringwidth > r)=1;

%% tile rings
tilenum=3;
ringimgs =[];
for k=1:tilenum
    ringimgs=[ringimgs ringimg];
end

ringimg=ringimgs;
ringimgs=[];

for k=1:tilenum
    ringimgs=[ringimgs;ringimg];
end

% imagesc(ringimgs);
% colormap bone;

filled_fullneon=zeros(screen_height_in_pixels,screen_width_in_pixels);
filled_fullneon(new_stimrect(2)+1:new_stimrect(4),new_stimrect(1)+1:new_stimrect(3))=ringimgs;
% imagesc(fullneon)
[x,y]=meshgrid(-ringsize:ringsize,-ringsize:ringsize);
[~,r]=cart2pol(x,y);
ringimg=x*0+.5;

rings=neon_exp.rings;
for ring=rings
    ringimg(ring-ringwidth < r & ring+ringwidth > r)=1;
end


%% tile rings
tilenum=3;
ringimgs =[];
for k=1:tilenum
    ringimgs=[ringimgs ringimg];
end

ringimg=ringimgs;
ringimgs=[];

for k=1:tilenum
    ringimgs=[ringimgs;ringimg];
end

fullneon=zeros(screen_height_in_pixels,screen_width_in_pixels);
fullneon(new_stimrect(2)+1:new_stimrect(4),new_stimrect(1)+1:new_stimrect(3))=ringimgs;

