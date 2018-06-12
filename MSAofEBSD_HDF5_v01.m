%% Matlab script for factor/cluster analysis of EBSD data
%                                  - Angus J Wilkinson May 2018

% clear the decks by removing all variables and shutting all matlab figures
clear variables; close all;

%start a timer - so we can tell how long different tasks took
tic;

%%  <<<< USER INPUTS >>>>
%   give location of dataset and specify pattern size to work on 
%   array size can get very large so may need to resize to bin down from
%   experimental conditions


% small ferritic steel sample - from Merlin- from DMC
BaseFolder='D:\Users\AJW\EBSD\Steel\Merlin\';
HDF5_filename='DX54_biaxial_ref_small.h5';
seg_method=3; %0-PCA, 1-PCA+Varimax, 3-kmeans
% these are user input - recorded patterns will be binned down to this size
PatSize2=800/8;  % generally need to reduce this from as recorded
PatSize1=576/8;  % generally need to reduce this from as recorded
NumPCA=22;  % user input of number of basis patterns
% select method for removing background 
flatfield=1; sigma=PatSize1/10.;  % Flatfield = 0-no removal, 1-divide by gausian
% blurred bckgnd for each image, 2-divide by mean over all images, 3-divide
% by recorded static background 
% 


%%  Setting up to read in data

% setup folder to save results in
if seg_method==0
    SaveFolder=[BaseFolder 'PCAfromHDF5'];
elseif seg_method==1
    SaveFolder=[BaseFolder 'VMfromHDF5'];
elseif seg_method==3
    SaveFolder=[BaseFolder 'KMfromHDF5'];
end    
mkdir(SaveFolder);
SaveFolder=[SaveFolder '\'];
mkdir([SaveFolder 'Patterns']);
mkdir([SaveFolder 'Maps']);

% open HDF5 file and read in some data on EBSD map
InputUser.HDF5_folder=BaseFolder;
InputUser.HDF5_file=HDF5_filename;
[ MapData,MicroscopeData,PhaseData,EBSP ]=bReadHDF5( InputUser );

% this bit finds the extent of the cropped region of the SEM image that was
% actually mapped with EBSD
MapSize1=max(MapData.YBeam)-min(MapData.YBeam)+1; % number of rows
MapSize2=max(MapData.XBeam)-min(MapData.XBeam)+1; % number of columns
MicroscopeData.NPoints;
Image=transpose(MicroscopeData.SEMImage(:,:,1));
figure;imagesc(Image);axis image off;colormap('gray');
figure;imagesc(transpose(reshape(MapData.RadonQuality,MapSize2,MapSize1)));axis image off;colormap('gray');

elapsedTime=toc;
hours=floor(elapsedTime/(60*60));
mins=floor((elapsedTime-60*60*hours)/(60));
secs=elapsedTime-60*60*hours - 60*mins;
display(['now setup to start reading data, time: hours=' num2str(hours) '  mins=' num2str(mins) '  secs=' num2str(secs)]);

    
%%  setup PCA analysis of vectorised patterns - select only strong contrast region     
% data is in a sequence of patterns, with associated sequence of positions, Euler angles etc 
testArray=zeros(PatSize1*PatSize2,MicroscopeData.NPoints);

for pattern_number=1:MicroscopeData.NPoints  %MapSize1*MapSize2
    RefPat=(double(h5read(EBSP.HDF5_loc,EBSP.PatternFile,double([1 1 pattern_number]),double([EBSP.PW EBSP.PH 1])))');
    RefPat=medfilt2(RefPat,'symmetric');
%     RefPat=imresize(RefPat,[PatSize1, PatSize2]); % scale the image size down
if flatfield==1
    BckGnd=  imgaussfilt(RefPat,sigma*size(RefPat,1)/PatSize1);
%    RefPat=medfilt2(RefPat,'symmetric');
    RefPat=RefPat./BckGnd;
    BckGnd=  imgaussfilt(RefPat,sigma*size(RefPat,1)/PatSize1);
    RefPat=RefPat-BckGnd;
    RefPat=medfilt2(RefPat,'symmetric');
end
RefPat=imresize(RefPat,[PatSize1, PatSize2]); % scale the image size down
testArray(:,pattern_number)=reshape(RefPat,PatSize1*PatSize2,1);
    if ( pattern_number==1000*floor(pattern_number/1000) )
        display(['  got to pattern number ',num2str(pattern_number), '   out of ',num2str(MicroscopeData.NPoints)]);
    end
    row(pattern_number)=MapData.YBeam(pattern_number)-min(MapData.YBeam)+1; %Y position in cropped region actually mapped with EBSD
    col(pattern_number)=MapData.XBeam(pattern_number)-min(MapData.XBeam)+1; %X position in cropped region actually mapped with EBSD
end

elapsedTime = toc;
hours=floor(elapsedTime/(60*60));
mins=floor((elapsedTime-60*60*hours)/(60));
secs=elapsedTime-60*60*hours - 60*mins;
display(['loaded all patterns into testarray, time: hours=' num2str(hours) '  mins=' num2str(mins) '  secs=' num2str(secs)]);

if flatfield==2
    % calc mean intensity - plus gaussian filter to smooth it
    meanPat=mean(testArray,2);
    meanPat=reshape(mean(testArray,2),PatSize1,PatSize2);
    meanPat= imgaussfilt(meanPat, 3);
    meanPat=reshape(meanPat,PatSize1*PatSize2,1);
end
if flatfield==3
    % read in static background, resize and divide each test pattern
    BckPat=(double(h5read('W:\AJW\PCA tests\Merlin\DX54_biaxial_ref_small.H5','/Scan 0/EBSD/Header/StaticBackground',double([1 1 1]),double([320 240 1])))');
    BckPat=imresize(BckPat,[PatSize1, PatSize2]); % scale the image size down
    meanPat=reshape(BckPat,PatSize1*PatSize2,1);
end
if (flatfield==2 || flatfield==3)
    for i=1:MicroscopeData.NPoints
        testArray(:,i)=testArray(:,i)./meanPat; % divide by background
%         testArray(:,i)=testArray(:,i)-meanPat; % subtract by background
    end
end
%%  run PCA analysis to split data into basis patterns (eigenvectors) and signal strength in map (loadings)
%   this uses a user input choice on number of factors to retain (NumPCA)
if seg_method<=1
    [coeff,score,latent,~,explained]=pca(testArray, 'Centered',false, 'NumComponents',NumPCA);
    % [coeff,score,latent,~,explained]=pca(testArray,'NumComponents',NumPCA);
    
    elapsedTime=toc;
    hours=floor(elapsedTime/(60*60));
    mins=floor((elapsedTime-60*60*hours)/(60));
    secs=elapsedTime-60*60*hours - 60*mins;
    display(['first pass of PCA on data, time: hours=' num2str(hours) '  mins=' num2str(mins) '  secs=' num2str(secs)]);
    
    grainTotal=zeros(MapSize1,MapSize2);
    grainMax=zeros(MapSize1,MapSize2);
    grainNext=zeros(MapSize1,MapSize2);
    grainIndex=zeros(MapSize1,MapSize2);
    for j=1:NumPCA
        pcaPat{j}=reshape(score(:,j),PatSize1,PatSize2);
        %    grain50{j}=reshape(coeff(:,j),MapSize2,MapSize1)';
        for k=1:MicroscopeData.NPoints
            grain50{j}(row(k),col(k))=coeff(k,j);
        end
        grainTotal=grainTotal+abs(grain50{j});
        tag=j*ones(MapSize1,MapSize2);
        grainIndex(abs(grain50{j})>grainMax)=tag(abs(grain50{j})>grainMax);
        grainNext(abs(grain50{j})>grainMax)=abs(grainMax(abs(grain50{j})>grainMax));
        grainMax(abs(grain50{j})>grainMax)=abs(grain50{j}(abs(grain50{j})>grainMax));
        
    end
    
    if seg_method==1
        % rotate basis vectors to find the VARIMAX solution
        [coeffVM, RotVM] = rotatefactors(coeff(:,1:NumPCA),'Method','varimax', 'Maxit', 5000, 'Reltol', 0.5*sqrt(eps));
        %[coeffVM, RotVM] = rotatefactors(coeff(:,1:NumPCA),'Method','varimax');
        elapsedTime=toc;
        hours=floor(elapsedTime/(60*60));
        mins=floor((elapsedTime-60*60*hours)/(60));
        secs=elapsedTime-60*60*hours - 60*mins;
        display(['hyper-rotate factors to VARIMAX solution, time: hours=' num2str(hours) '  mins=' num2str(mins) '  secs=' num2str(secs)]);
    elseif seg_method==0
        coeffVM=coeff(:,1:NumPCA);
        RotVM=eye(NumPCA);
    end
    scoreVM=score(:,1:NumPCA)*RotVM;
    grainVMtotal=zeros(MapSize1,MapSize2);
    grainVMmax=zeros(MapSize1,MapSize2);
    grainVMnext=zeros(MapSize1,MapSize2);
    grainVMindex=zeros(MapSize1,MapSize2);
    grainVMsign=zeros(MapSize1,MapSize2);
    
    % reform the basis patterns, and make maps showing their signal strength
    for j=1:NumPCA
        pcaVMpat{j}=reshape(scoreVM(:,j),PatSize1,PatSize2);
        %     grainVM{j}=reshape(coeffVM(:,j),MapSize2,MapSize1)';
        for k=1:MicroscopeData.NPoints
            grainVM{j}(row(k),col(k))=coeffVM(k,j);
        end
        grainVMtotal=grainVMtotal+abs(grainVM{j});
        tag=j*ones(MapSize1,MapSize2);
        grainVMindex(abs(grainVM{j})>grainVMmax)=tag(abs(grainVM{j})>grainVMmax);
        grainVMnext(abs(grainVM{j})>grainVMmax)=abs(grainVMmax(abs(grainVM{j})>grainVMmax));
        grainVMmax(abs(grainVM{j})>grainVMmax)=abs(grainVM{j}(abs(grainVM{j})>grainVMmax));
    end
    
    % check for negative contrast and invert relevant patterns
    % save basis patterns as tif files
    for j=1:NumPCA
        signPat(j)=sign(mean(grainVM{j}(grainVMindex==j)));
        meanPatPCA(j)=mean2(grain50{j});
        meanPatVM(j)=mean2(grainVM{j});
        pcaVMpat{j}=pcaVMpat{j}*signPat(j);
        grainVM{j}=grainVM{j}*signPat(j);
        %     MinInt=min(scoreVM(:,j));
        %     MaxInt=max(scoreVM(:,j));
        
        MinInt2=min(min(grainVM{j}));
        MaxInt2=max(max(grainVM{j}));
        imwrite(((grainVM{j}-MinInt2))/(MaxInt2-MinInt2),[SaveFolder 'Maps\map_' num2str(j,'%5.5u') '.tif'])
        
        MinInt=min(min(pcaVMpat{j}));
        MaxInt=max(max(pcaVMpat{j}));
        imwrite(((pcaVMpat{j}-MinInt))/(MaxInt-MinInt),[SaveFolder 'Patterns\pattern_' num2str(j,'%5.5u') '.tif'])
    end
    
    figure(303);imagesc(grainVMtotal);axis image;colormap('jet');colorbar;title('Total Signal');saveas(303,[SaveFolder 'VM total signal.tif']);
    figure(304);imagesc(grainVMmax);axis image;colormap('jet');colorbar;title('Strongest Signal');saveas(304,[SaveFolder 'VM strongest signal.tif']);
    figure(305);
    subplot(1,2,1);imagesc(grainVMindex);axis image off;colormap('jet');
    subplot(1,2,2);histogram(grainVMindex(isfinite(grainVMindex)),0:NumPCA);saveas(305,[SaveFolder 'VM Grain Indexing.tif']);
    figure(306);imagesc(grainVMmax./grainVMtotal);axis image off;colormap('jet');colorbar;title('ratio strongest to total signal');saveas(306,[SaveFolder 'VM ratio strongest to total signal.tif']);
    figure(307);imagesc(grainVMnext);axis image off;colormap('jet');colorbar;title('second strongest signal');saveas(307,[SaveFolder 'VM second strongest signal.tif']);
    figure(308);imagesc(log10(grainVMmax./grainVMnext));axis image off;colormap('jet');colorbar;title('log ratio of first and second strongest signal');caxis([0 3]);saveas(308,[SaveFolder 'VM log ratio of first and second strongest signal.tif']);
    figure(309);imagesc(grainVMindex);axis image off;colormap('jet');colorbar;title('segmented grains');saveas(309,[SaveFolder 'VM segmented grains.tif']);
    
    
    elapsedTime = toc;
    hours=floor(elapsedTime/(60*60));
    mins=floor((elapsedTime-60*60*hours)/(60));
    secs=elapsedTime-60*60*hours - 60*mins;
    display(['all done, time: hours=' num2str(hours) '  mins=' num2str(mins) '  secs=' num2str(secs)]);
    
    %
    first_few=20;
    if NumPCA<=first_few
        first_few=NumPCA;
    end
    for j=1:first_few
        figure(401);subplot(4,5,j);imagesc(pcaVMpat{j});axis image off;colormap('gray');
        figure(402);subplot(4,5,j);imagesc(grainVM{j});axis image off;colormap('gray');
    end
    saveas(401,[SaveFolder 'VM Patterns inversion corrected.tif']);
    saveas(402,[SaveFolder 'VM Grains inversion corrected.tif']);
    
end

%% use k-means clustering to segment EBSD data
if seg_method==3

    [idx,KMcentroid,sumd, Dist] = kmeans(testArray',NumPCA,'Replicates',50); 
    KMgrainMap=transpose(reshape(idx,MapSize2,MapSize1));
    [Dlist,KMindex]=sort(Dist,2);
    KMratioMap=transpose(reshape((Dlist(:,2)./Dlist(:,1)),MapSize2,MapSize1));
    
    elapsedTime = toc;
    hours=floor(elapsedTime/(60*60));
    mins=floor((elapsedTime-60*60*hours)/(60));
    secs=elapsedTime-60*60*hours - 60*mins;
    display(['k-means cluster analysis done: hours=' num2str(hours) '  mins=' num2str(mins) '  secs=' num2str(secs)]);
    
    figure(608);imagesc(KMratioMap);axis image off;colormap('jet');colorbar; title('k-means grain ratio of distances');
    saveas(608,[SaveFolder 'KM ratio of distances.tif']);
    
    figure(609);imagesc(KMgrainMap);axis image off;colormap('jet');colorbar; title('k-means grain segmentation');
    saveas(609,[SaveFolder 'KM segmented grains.tif']);

    
    MaxDist=max(max(Dist));
    MinDist=min(min(Dist));
    for j=1:NumPCA
        KMpat{j}=reshape(KMcentroid(j,:),PatSize1,PatSize2);
        KMdist{j}=transpose(reshape(Dist(:,j),MapSize2,MapSize1));
        if j<=20
            figure(601);subplot(4,5,j);imagesc(KMpat{j});axis image off;colormap('gray');
            figure(602);subplot(4,5,j);imagesc(MaxDist-KMdist{j});axis image off;colormap('gray');caxis([0 MaxDist-MinDist]);
        end
%         MinInt2=min(min(KMdist{j}));
%         MaxInt2=max(max(KMdist{j}));
        MaxInt2=MaxDist;MinInt2=MinDist;
        imwrite(((MaxInt2-KMdist{j}))/(MaxInt2-MinInt2),[SaveFolder 'Maps\map_' num2str(j,'%5.5u') '.tif'])
        
        MinInt=min(min(KMpat{j}));
        MaxInt=max(max(KMpat{j}));
        imwrite(((KMpat{j}-MinInt))/(MaxInt-MinInt),[SaveFolder 'Patterns\pattern_' num2str(j,'%5.5u') '.tif'])

    end
    saveas(601,[SaveFolder 'KM Patterns inversion corrected.tif']);
    saveas(602,[SaveFolder 'KM Grains inversion corrected.tif']);

%     figure;
%     silhouette(testArray',idx);
%     E = evalclusters(testArray','kmeans','silhouette','klist',[10:20]);
%     
end


%% save some output files

if seg_method==3
    grainMap=KMgrainMap;
elseif seg_method==1
    grainMap=grainVMindex;
elseif seg_method==0
    grainMap=grainVMindex;
end    
save([SaveFolder 'GrainMap.mat'],'MicroscopeData', 'grainMap' );

