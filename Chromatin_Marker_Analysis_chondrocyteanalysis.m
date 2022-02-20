function Chromatin_Marker_Analysis_FLEX_v2
clc;  close all; clear all;
%% reads in immunostaining data with different channels for spatial analysis

filename_avg=['Compiled_data']; %% Change everytime! 
stored_path=evalin('base','path'); % get current data filepath (if stored before)
home_path=cd; % get path of matlab file
[im_name,im_path]=uigetfile('*.nd2','Choose Files','MultiSelect','on',stored_path); % load names and pathes of files
sample_num=length(im_name); % number of samples in the sweep stack
filename_save=[im_path 'MarkerAnalysis2.xlsx'];

 
%% Structure Initiations
nuc_props_save = struct;
dist_analysis= struct;
k9_foci=struct;
foci_struct=struct;

%% LOOP

for z=1:length(im_name)

% ANALYSIS OPTIONS
distance_analysis=1;
    bin_num=5;
    bin_num_2=3;
    bin_num_3=10;

    threshold_method='hist'; % choose threshold method for segmentation: 'hist' (using defined cut of in histogram) , 'gray' (graythresh)
    hist_thresh=66;   
    fig=1;

%% SETTINGS
file_format='nd2'; % choose between 8-bit "tif" or "nd2" as image type (nd2 is prefered)

bin_size=3; %bin size for chromatin histogram, default=3

binary_thresh_factor=[1 1 1 1]; % modifies graythresh-derived threshold value for binary segmentation of each color channel, default=1
OEC=[-1 -1 -1 -1]; % over-exposure correction in % (gets additionally substracted from 99%)

smooth_border=1; % smooth nuclear border dilation and erosion; 1=on
dilate_border=0; % extend border outline, default=0
threshold_enhancer=.1; % factor by which threshold value for border detection is reduced; 1=no change

%colormap:
H=hsv(100); J=jet(100); cmp=[0 0 0; 1 1 1; J(30,:); H(30,:); J(90,:); H(88,:)];
%% IMAGE OPTIONS

ref_channel=1; % used to find nuclear boundry
marker_channel= [1,3,4]; % channels to be analyzed
marker_name={'DAPI','H3K27','H3K9'}; % choose names for each channel to display in excel afer saving
H3K9Channel=3; %CHANGE IF H3K9 channel is stored differently!! H3K19 # in MC

%% INPUT

filename=strcat(im_name{z});
data=bfopen([im_path filename]); 
D=data{1,1}; % convert data from nd2 file to tif
assignin('base','path',path);
stored_path=evalin('base','path');

% read in reference channel
for r=1:length(ref_channel);
    RC(:,:,r)=D{ref_channel(r),1};
end

mc_num=length(marker_channel); % number of marker channels
for c=1:mc_num
    MC(:,:,c)=D{marker_channel(c),1}; % marker channels
end

omeMeta = data{1, 4};    
% returns the default unit type
voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm/pixel in x
resf=double(voxelSizeX);
assignin('base','path',path); % store current filepath for next execution in workspace

%% Segment Nucleus

[segmented_image, final_mask_image]=Image_Segment(imadjust(RC(:,:,1)));
IM_processed=segmented_image;
rc_bw=final_mask_image;
rc_bw=bwareaopen(rc_bw,100);

[bound,label]=bwboundaries(rc_bw); % get nuclear boundaries
nuc_num=length(bound); % number of identified nuclei

if nuc_num<1;
    rc_bw=im2bw((RC(:,:,1)),0.01);
    rc_bw=imfill(rc_bw,'holes'); % fill holes
    rc_bw=bwareaopen(rc_bw,100); % filter small objects
    [bound,label]=bwboundaries(rc_bw); % get nuclear boundaries
    nuc_num=length(bound);    % number of identified nuclei
end

prop=regionprops(rc_bw,'Centroid','Area','MajorAxisLength','Perimeter','BoundingBox'); % get region properties


% OVERLAY BOUNDARY AND NUMBER NUCLEI IN ORIGINAL IMAGE
figure; imshow(RC(:,:,1)); caxis auto; hold on;
for p=1:nuc_num
    border=bound{p};
    cent(p,:)=prop(p).Centroid; % save into seperate centroid varibale
    area(p)=prop(p).Area;
    peri(p)=prop(p).Perimeter;
    plot(border(:,2),border(:,1),'r')
    text(cent(p,1),cent(p,2),num2str(p),'color','r');
end

% FILTER AREAs
circ=4*pi*area./(peri.^2); % circularity of areas
del=find(circ<.2); % index for areas to be deleted based on low circularity

n=find(circ>.2);

for d=1:length(del)
    border=bound{del(d)};
    plot(border(:,2),border(:,1),'b')
end

%choose center most nucleus out of the most circular nuclei  
if length(n)>1
for cirnuc=1:length(n);
    cir_cent(cirnuc,:)=prop(n(cirnuc)).Centroid; 
end
%chose most center nucleus
dxy = fliplr(size(rc_bw)/2) - vertcat(cir_cent);  %distance along x and y from centre to centroid
[~, mostcentral] = min(hypot(dxy(:, 1), dxy(:, 2)));  %most central is the index of the object the closest to the centre
n=n(mostcentral);

end

bw=rc_bw;

%Nuclear Properties
%Find nuclear actual area in um^2
prop=regionprops(bw,'Centroid','PixelList','Orientation','Area','MajorAxisLength','MinorAxisLength','Perimeter');
area_umsq=(prop(n).Area)*(double(voxelSizeX))^2;
majoraxis=prop(n).MajorAxisLength*double(voxelSizeX);
minoraxis=prop(n).MinorAxisLength*double(voxelSizeX);
aspectratio=majoraxis/minoraxis;
nuclear_perimeter=prop(n).Perimeter;
%Save nuclear properties to Structure when in loop

nuc_props_save(z).name = filename;
nuc_props_save(z).nuclear_area = area_umsq;
nuc_props_save(z).nuclear_aspratio = aspectratio;

% update all varibale to choosen nucleus
border=bound{n};
nuc_props=prop(n);
c=prop(n).Centroid; nuc_cent(1)=c(2); nuc_cent(2)=c(1);

% get all pixel coordinates of selected nucleus
temp=nuc_props.PixelList; % pixel list of current nucleus, save in temporary varibale
nuc_pix=[temp(:,2), temp(:,1)]; % all coordinates of pixels in the nucleus

bw_nuc=zeros(size(MC(:,:,1)));
for i=1:length(nuc_pix)
bw_nuc(nuc_pix(i,1),nuc_pix(i,2))=1;
end

%% (2) CHANNEL SEGMENTATIONS

% Generate Nucleus-Only Maps for each Marker Channel

%convert image informtion into linear vectors
for m=1:mc_num
    for i=1:length(nuc_pix)
        LINEAR(i,m)=double(MC(nuc_pix(i,1),nuc_pix(i,2),m));
    end
    % pre-normalize linear color vectors [0 1]
    LIN_norm(:,m)=LINEAR(:,m)-min(LINEAR(:,m)); % baseline substraction
    LIN_norm(:,m)=LIN_norm(:,m)./max(LIN_norm(:,m)); % normalize
    
    % histogram normalize to avoid from few high intense pixels
    % --> plateau highest 1% of pixels
    
    [count(:,m),center(:,m)]=hist(LIN_norm(:,m),100); % generate intensity histograms
    % convert to summed histogram
    hsum(:,m)=count(:,m);
    for i=2:length(count(:,m)); hsum(i,m)=hsum(i-1,m)+hsum(i,m); end
    hsum(:,m)=hsum(:,m)./max(hsum(:,m)).*100; % normalize to [0 100%]
    
    % find 99% cut-offs in summed histogram, OEC=over-exposure correction
    [idx_temp]=knnsearch(hsum(:,m),99-OEC(m)); CO99(m)=center(idx_temp,m);
    % find cut off for thresholding
    [idx_temp]=knnsearch(hsum(:,m),hist_thresh); COth(m)=center(idx_temp,m); % cut-off threshold
    
    % cut of 1% of the intensity peak through renormalization to 99% cut-off
    LIN_norm(:,m)=LIN_norm(:,m)./CO99(m); LIN_norm(LIN_norm(:,m)>1,m)=1;
end

figure(fig); fig=fig+1;
% plot histogram normalization
for m=1:mc_num
    subplot(mc_num,3,(m-1)*3+1:(m-1)*3+2);
    plot(center(:,m),hsum(:,m),'color',cmp(marker_channel(m)+2,:),'linewidth',2); hold on; line([COth(m),COth(m)],[0 100],'Color','k')
    % subplot for adjusted histograms see below
end


%% BINARY SEGMENTATION
% Determine pixels that are positive or negative for the respective color
% channel. One pixel can be positive for multiple channels. Segregation is
% achieved through binary thresholding using the GRAYTHRESH function plus a
% factor, binary_thresh_factor, for manual adjustment; default=2
% 
% generate nuclear-only contineous (IM_nuc) and binary (BW_nuc) maps from normalized line vectors
% Also generate new line vectors for binary maps
IM_nuc=zeros(size(MC)); % normalized nucleus only image matrix
BW_nuc=zeros(size(MC)); % binary nucleus only map
BW_lin=zeros(length(nuc_pix),mc_num); % binary nucleus only vector
void=zeros(size(MC(:,:,1))); % map for "empty" pixels that have not been assigned to any category
void_lin=zeros(length(nuc_pix),1);
% 
for m=1:mc_num
    for n=1:length(nuc_pix)
        IM_nuc(nuc_pix(n,1),nuc_pix(n,2),m)=LIN_norm(n,m);
    end
    % set threshold for segmentation
    if strcmp(threshold_method,'gray')==true
        temp=IM_nuc(:,:,m); temp(temp==0)=[]; bin_thresh(m)=graythresh(temp); % load current nuc only map into temporary variable and delete all 0 entries for proper threshold calculation
    end
    if strcmp(threshold_method,'hist')==true
        bin_thresh(m)=COth(m);
    end
    % generate binary
    BW_nuc(:,:,m)=im2bw(IM_nuc(:,:,m),bin_thresh(m)*binary_thresh_factor(m));
    
    % Plot adjusted histograms with threshold lines
    fig=fig-1; % draw threshold lines into previous figure
    figure(fig); fig=fig+1;
    subplot(mc_num,3,m*3)
    h1=histogram(LIN_norm(:,m));
    line([bin_thresh(m),bin_thresh(m)],[0 max(h1.Values)],'color','k','linewidth',2); h1.FaceColor=cmp(marker_channel(m)+2,:);
    
    % plot adjusted and binary nuclear maps of all marker channels
    figure(fig);
    subplot(2,mc_num,m); title(['Channel ' num2str(m)])
    imshow(IM_nuc(:,:,m)); hold on; plot(border(:,2),border(:,1),'r')
    subplot(2,mc_num,mc_num+m)
    imshow(BW_nuc(:,:,m)); hold on; plot(border(:,2),border(:,1),'r')
end
% fig=fig+1;

    msc_count=zeros(4); % initialize marker single coverage count variable
    for i=1:length(nuc_pix)
        
        void_check=true; % checks if current location is empty
        for m=1:mc_num
            if BW_nuc(nuc_pix(i,1),nuc_pix(i,2),m)==1; BW_lin(i,m)=1; void_check=false; end
        end
        if void_check==true; void(nuc_pix(i,1),nuc_pix(i,2))=1; void_lin(i)=1; else
            % collect number of points for marker single coverage (MSC)
            if sum(BW_lin(i,:))==1 % if exactly one channel is present
                idx=find(BW_lin(i,:)==1); % find which channel
                msc_count(idx)=msc_count(idx)+1; % increase its single count
            end
        end
    end
%     
    % VISUAL OUTPUT
    figure(fig); fig=fig+1; imshow(void); hold on; % show composite map including "void" pixels
    % overlay outlines of other areas
    for m=1:mc_num
        bound=bwboundaries(BW_nuc(:,:,m));
        for i=1:length(bound); b=bound{i}; plot(b(:,2),b(:,1),'color',cmp(marker_channel(m)+2,:),'linewidth',2); end
    end   

        
%% Calculate Foci H3K9 channel
list(z).name = filename;
resize_factor=1; % factor for image resizing before peak analysis
    focis=IM_nuc(:,:,H3K9Channel);
    focis=imresize(focis,resize_factor); % resize image for better peak detection
    focis_bw=imresize(bw_nuc,resize_factor);%im2bw(focis,0.01);
    peri_b=bwboundaries(focis_bw,'noholes');
    stats_foci = regionprops(focis_bw,'Centroid','MinorAxisLength');
    [pk] = FastPeakFind(focis); % find peaks
    pkx3=pk(1:2:end); pky3=pk(2:2:end);
    %foci3=length(pkx3); % number of  foci
    figure,imshow(focis); hold on 
    plot(pk(1:2:end),pk(2:2:end),'r+')
    for p=1:length(peri_b);
        bound_p=peri_b{p};
        plot(bound_p(:,2),bound_p(:,1),'b')
    end
    perimeter_yx=peri_b{1,1};
    perimeter_xy=[perimeter_yx(:,2),perimeter_yx(:,1)];
    x_foci=pk(1:2:end);
    y_foci=pk(2:2:end);
    foci_center_count=0;

    for f_c=1:length(x_foci)

        pto1 = [stats_foci.Centroid(1) stats_foci.Centroid(2)];  
        pto2 = [x_foci(f_c) y_foci(f_c)];
        V = pto2 - pto1;
        % The distance between the points would be:
        distance_from_centroid(f_c) = norm(V);
        % which will be extended (by 20% in this case) here
        factor_distance = 500;
        % Extend the ray
        pext = pto1 + V*factor_distance;
        % plot...
        %plot([pto1(1),pto2(1)],[pto1(2),pto2(2)],'bo',[pto1(1),pext(1)],[pto1(2),pext(2)],'r-')
        line([x_foci(f_c),stats_foci.Centroid(1)],[y_foci(f_c),stats_foci.Centroid(2)],'LineWidth',1);
        %hLine=line([x_foci(f_c),pext(1)],[y_foci(f_c),pext(2)],'Color','red');

        hLine = imline(gca,[x_foci(f_c),pext(1)],[y_foci(f_c),pext(2)]); %[x1 y1; x2 y2].
        singleLineBinaryImage = hLine.createMask();
        se = strel('square',2);
        singleLineBinaryImage = imdilate(singleLineBinaryImage,se);
        [y_line, x_line] = find(singleLineBinaryImage);
        line_points=[x_line, y_line];
        [logical, loca_in_peri] = ismember(line_points,perimeter_xy,'rows');
        indx_intersect=find(logical);
        point_intersect=line_points(indx_intersect,:);

        plot(point_intersect(1,1),point_intersect(1,2),'g*');
        pto3=[point_intersect(1,1) point_intersect(1,2)];
        V2=pto3-pto2;
        distance_from_foci = norm(V2);
        distance_ratio(f_c) = distance_from_centroid(f_c)/(distance_from_centroid(f_c)+distance_from_foci);

        if distance_from_centroid(f_c)<((stats_foci.MinorAxisLength/2)-(.1*stats_foci.MinorAxisLength));
            foci_center_count=foci_center_count+1;
        end
    end

    number_foci=foci_center_count;
    foci_struct_FPF(z).K9number=number_foci;
    avg_dist_ratio(z)=mean(distance_ratio);
    dist_ratio(z).all=distance_ratio;
            

%% (3) MARKER ARCHITECUTRE ANALYSIS

% info: coordinates 1=Y, 2=X

% get distance along border from nuc center in relation to its angle phi [0 360]
% --> Use to normalize every pixel to total length based on its angel
for b=1:length(border)
    dX=border(b,2)-nuc_cent(2); dY=-1*(border(b,1)-nuc_cent(1));
    border_distance(b)=sqrt(dX^2+dY^2);
    phi(b)=atan(dY/dX);

    % adjust angle to [1 360]
    if dX<=0; phi(b)=phi(b)+pi; end
    if dY<=0 && dX>0; phi(b)=phi(b)+2*pi; end
end

phi=phi./pi.*180; % convert to rad
% remove redundancy and sort
[phi,I]=unique(phi); border_distance=border_distance(I);
% interpolation    
rad=1:360; % angle-vector for interpolation
border_distance_int=spline(phi,border_distance,rad); % interpolated border distance for easier recall

% calculate relative distance and angle beta [0 90] between minor and major axis
for p=1:length(nuc_pix)
    dX=nuc_pix(p,2)-nuc_cent(2); dY=-1*(nuc_pix(p,1)-nuc_cent(1));
    pix_distance=sqrt(dX^2+dY^2); % distance from nuc center [pix]
    beta90(p)=atan(dY/dX); % angle from nuc center [RAD]

    % adjust angle to [1 360]
    if dX<=0; beta360(p)=beta90(p)+pi; else
        if dY<=0 && dX>0; beta360(p)=beta90(p)+2*pi; else
            beta360(p)=beta90(p);
        end
    end
    beta360(p)=beta360(p)./pi.*180; beta360(p)=round(beta360(p)); % convert to DEG from RAD and round
    if beta360(p)==0; beta360(p)=360; end

    % normalize pixel distance using angle-distance relationship from above
    pix_dist_rel(p)=pix_distance./border_distance_int(beta360(p)); % relative distance of pixel from nuc center [0 1]
end
beta90=round(beta90./pi.*180); % convert to DEG from RAD and round

% DISTANCE ANALYSIS
if distance_analysis==1
    
    % SEGMENTED
    % split distance and angle data for each channel seperatly
    % NOTE: Channel overlapp depends on settings on top!
    xbins_100=0.005:0.01:1;
    xbins_n= 0.5/bin_num:1/bin_num:1;
    xbins_n2=0.5/bin_num_2:1/bin_num_2:1;
    xbins_n3=0.5/bin_num_3:1/bin_num_3:1; 
    
    for m=1:mc_num
        DIST(m)={pix_dist_rel(BW_lin(:,m)==1)};
        DIST_100(:,m)=hist(DIST{m},xbins_100);
        DIST_n(:,m)=hist(DIST{m},xbins_n);
        DIST_n2(:,m)=hist(DIST{m},xbins_n2);
        DIST_n3(:,m)=hist(DIST{m},xbins_n3);
    end
    % void channel
    DIST(mc_num+1)={pix_dist_rel(void_lin==1)};
    DIST_100(:,mc_num+1)=hist(DIST{mc_num+1},xbins_100);
    DIST_n(:,mc_num+1)=hist(DIST{mc_num+1},xbins_n);
    DIST_n2(:,mc_num+1)=hist(DIST{mc_num+1},xbins_n2);
    DIST_n3(:,mc_num+1)=hist(DIST{mc_num+1},xbins_n3);
    

    % convert distance to relative values (% of distance bin)
    for d=1:length(xbins_100); DIST_100(d,:)=DIST_100(d,:)./sum(DIST_100(d,:)); end
    for d=1:length(xbins_n); DIST_n(d,:)=DIST_n(d,:)./sum(DIST_n(d,:)); end
    for d=1:length(xbins_n2); DIST_n2(d,:)=DIST_n2(d,:)./sum(DIST_n2(d,:)); end
    for d=1:length(xbins_n3); DIST_n3(d,:)=DIST_n3(d,:)./sum(DIST_n3(d,:)); end
   
    figure(fig); fig=fig+1; 
    subplot(1,2,1); hold on;
    for m=1:mc_num
        plot(xbins_100,DIST_100(:,m),'color',cmp(marker_channel(m)+2,:),'linewidth',3)
    end
    plot(xbins_100,DIST_100(:,mc_num+1),'color',cmp(2,:),'linewidth',3)
    set(gca,'Color','k')
    xlabel('rel. nuclear center distance'); ylabel('vol. fraction')
    
    subplot(1,2,2)
    b=bar(DIST_n);
    for m=1:mc_num
        b(m).FaceColor=cmp(marker_channel(m)+2,:);
    end
    b(mc_num+1).FaceColor=cmp(2,:);
    set(gca,'Color','k')
    
    
end


% SAVE DATA to Excel
    % generate new filename
    filename(end-3:end)=[];
    filename_1=['\' filename 'Marker Analysis.xlsx'];

    cd(im_path)
    
    % DISTANCE
    if distance_analysis==1
        % write header
        header=cell(1,mc_num+2); % create header cell containing all headlines
        header(1)={'rel. nuc distance'}; header(2:end-1)=marker_name; header(mc_num+2)={'void'};
        xlswrite(filename_1,header,'distance','A1')
        header(1)={'n-bin:'};
        xlswrite(filename_1,header,'distance','G1')
        xlswrite(filename_1,header,'distance','M1')
        xlswrite(filename_1,header,'distance','S1')
        % write DISTANCE data
        RES=[xbins_100',DIST_100];
        xlswrite(filename_1,RES,'distance','A2')
        RES1=[xbins_n',DIST_n];
        xlswrite(filename_1,RES1,'distance','G2')
        RES2=[xbins_n2',DIST_n2];
        xlswrite(filename_1,RES2,'distance','M2')
        RES3=[xbins_n3',DIST_n3];
        xlswrite(filename_1,RES3,'distance','S2')
        NUC_SAVE=[area_umsq, aspectratio];
        header2={'nuclear area','nuclear aspect ratio'};
        xlswrite(filename_1,header2,'nuclear properties','A1')
        xlswrite(filename_1,NUC_SAVE, 'nuclear properties','A2')   
        
        dist_analysis(z).name = filename;
        dist_analysis(z).xbins_data = xbins_n';
        dist_analysis(z).distance_data = DIST_n;
        dist_analysis(z).xbins_data2 = xbins_n2';
        dist_analysis(z).distance_data2 = DIST_n2;
        dist_analysis(z).xbins_data3 = xbins_n3';
        dist_analysis(z).distance_data3 = DIST_n3;
        dist_analysis(z).xbins_data100 = xbins_100';
        dist_analysis(z).distance_data100 = DIST_100;
    end
    

close all
save('parameters','H3K9Channel', 'threshold_enhancer','im_name','im_path','nuc_props_save','dist_analysis','z','xbins_100','xbins_n','xbins_n2','xbins_n3','filename_avg','k9_foci','foci_struct','foci_struct_FPF','avg_dist_ratio','dist_ratio','list'); 
clear variables
load('parameters');

end

% Calculate average of parameters 
%convert structures to cells first
nuc_props_cell=struct2cell(nuc_props_save);
dist_analysis_cell=struct2cell(dist_analysis);
list_cell=struct2cell(list);
k9_foci_cell_FPR=struct2cell(foci_struct_FPF);

k9_foci_matFPR=cell2mat(k9_foci_cell_FPR(1,:,:));

nuc_area2=cell2mat(nuc_props_cell(2,:,:));
nuc_aspectratio2=cell2mat(nuc_props_cell(3,:,:));

FL=(list_cell(1,:,:)); 

distance_5=cell2mat(dist_analysis_cell(3,:,:));
distance_10=cell2mat(dist_analysis_cell(7,:,:));%% Figure out number 
distance_100=cell2mat(dist_analysis_cell(9,:,:));

avg_distance_5=mean(distance_5,3);
avg_distance_10=mean(distance_10,3);
avg_distance_100=mean(distance_100,3);
avg_nuclear_area=mean(nuc_area2,3);
avg_nuclear_aspectratio=mean(nuc_aspectratio2,3);

std_distance_5=std(distance_5,0,3);
std_distance_10=std(distance_10,0,3);
std_distance_100=std(distance_100,0,3);
std_nuclear_area=std(nuc_area2,0,3);
std_nuclear_aspectratio=std(nuc_aspectratio2,0,3);


%Save Compiled Data to Excel Workbook

header=cell(1,9); % create header cell containing all headlines
header(1)={'rel. nuc distance'}; header(2)={'DAPI'}; header(3)={'H3K27'};header(4)={'H3K9'};header(5)={'void'};header(6)={'Std DAPI'}; header(7)={'Std H3K27'};header(8)={'Std H3K9'};header(9)={'Std void'};
xlswrite(filename_avg,header,'distance','A1')
header(1)={'n-bin:'};
xlswrite(filename_avg,header,'distance','K1')
header(1)={'n-bin:'};
xlswrite(filename_avg,header,'distance','U1')

% write DISTANCE data
RES100=[xbins_100',avg_distance_100,std_distance_100];
xlswrite(filename_avg,RES100,'distance','A2')
RES5=[xbins_n',avg_distance_5,std_distance_5];
xlswrite(filename_avg,RES5,'distance','K2')
RES10=[xbins_n3',avg_distance_10,std_distance_10];
xlswrite(filename_avg,RES10,'distance','U2')

na=nuc_area2(1,1,:);
nar=nuc_aspectratio2(1,1,:);
filelist=string(FL(:));
NUC_AVGSAVE=[avg_nuclear_area,std_nuclear_area,avg_nuclear_aspectratio,std_nuclear_aspectratio];
header2={'nuclear area','std nuc area','nuclear aspect ratio','std nuc a.r.'};
xlswrite(filename_avg,header2,'nuclear properties','G1')
xlswrite(filename_avg,NUC_AVGSAVE, 'nuclear properties','G2')
Nuc_all=[filelist,na(:),nar(:)];
xlswrite(filename_avg,Nuc_all, 'nuclear properties','A1')

f_fpf=k9_foci_matFPR(1,1,:);
Foci_save_fpf=[filelist,f_fpf(:)];

header4={'filename','K9-Foci'};
xlswrite(filename_avg, header4, 'foci_fpf','C1');
xlswrite(filename_avg, Foci_save_fpf, 'foci_fpf','C2');
end
        
