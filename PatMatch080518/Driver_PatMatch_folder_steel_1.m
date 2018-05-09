%% This code creates a test library using MTex generated uniform odf
% library matching is then undertaken with all patterns in a folder
%

% clear existing
clear all;
close all;

       p = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(p)
            poolsize = 0;
        else
            poolsize = p.NumWorkers;
        end
        if poolsize <2
            parpool
        end

tic

%% Define Master Pattern and Detector Geometry
                  
crop_flag=0;
% Pattern=imread('Fe_2866_2866_ver3.tif'); % Import the Pattern File 
Pattern=imread('Fe_2866_2866_ver3.tif'); % Import the Pattern File 
Pattern=Pattern(:,:,1);Pattern = double(Pattern);
Pattern=imgaussfilt(Pattern,4);

sample_tilt=70.0;
camera_tilt=4.0;
tilt=(sample_tilt-90.0)-camera_tilt;

W = 100*4/4; % Width of the Detector Plane in Pixels
H = 72*4/4; % Height of the Detector Plane in Pixels



% steels patterns
Xstar=0.4771; % Pattern centre as a fraction of the Width from left of the Detector Plane[note Bruker use zero on left]
Ystar=0.7139; % Pattern centre as a fraction of the Height from bottom of the Detector Plane [note Bruker use zero at top]
D=0.8125*H/W; % Camera Length as a fraction of the Detector Width in Pixels [note Bruker use fraction of Height]


D = D*W/H;
points=500000;
InputFolder='D:\Users\AJW\EBSD\Steel\Merlin\KMfromHDF5_test22\';
GrainFile='GrainMap';
ExptPatFolder=[InputFolder 'Patterns\'];
OutputFolder=[InputFolder 'PMoutput\500k100_72rot\'];
mkdir(OutputFolder);


%% Create Structure

CS = crystalSymmetry('m-3m');
% SS = specimenSymmetry('orthorhombic');
SS = specimenSymmetry('triclinic');
% % plotting convention - this was AJW interpretation of 'which way is up' 
% setMTEXpref('xAxisDirection','east');
% setMTEXpref('zAxisDirection','outofPlane');

% plotting convention - this gives consistent plots with Bruker output
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');

% uniform crystal orientations
odf = uniformODF(CS,SS);

% % distribution near slected orientation - for testing
% RefOri = orientation('Euler',64.0058*degree,119.1247*degree,6.5392*degree,CS,SS); % Euler angles for steels pattern1
% 
%  psi = vonMisesFisherKernel('HALFWIDTH',20.0*degree);
%   odf = unimodalODF(RefOri,psi);


% calc random orientations from the odf distribution
ori = calcOrientations(odf,points);
count2=points;
% [phi1,Phi,phi2] = Euler(ori,'Bunge') ;
[phi1,Phi,phi2] = project2EulerFR(ori,CS,'Bunge');   
phi1 =(180.0/pi)*phi1;
Phi  = (180.0/pi)*Phi;
phi2 = (180.0/pi)*phi2;

%% Loop
    elapsedTime = toc;
    hours=floor(elapsedTime/(60*60));
    mins=floor((elapsedTime-60*60*hours)/(60));
    secs=elapsedTime-60*60*hours - 60*mins;
    display(['now setup to calculate pattern library: hours=' num2str(hours) '  mins=' num2str(mins) '  secs=' num2str(secs)]);

% Loop through each orientation
parfor i=1:count2
   Image = zxz_tilt_sim1([phi1(i),Phi(i),phi2(i)],[Xstar,Ystar,D],W,H,Pattern,tilt);
   Image = Convert_Image(Image,H,W);
   Image=uint8(Image);
   LibPattern{i}=Image;
end

% end
    elapsedTime = toc;
    hours=floor(elapsedTime/(60*60));
    mins=floor((elapsedTime-60*60*hours)/(60));
    secs=elapsedTime-60*60*hours - 60*mins;
    display(['pattern library calculated: hours=' num2str(hours) '  mins=' num2str(mins) '  secs=' num2str(secs)]);


%% Correllation between Experimental patterns and Simulated library
    
    FileList=dir(ExptPatFolder);
    FileNames=cell(size(FileList,1),1);
    FileAllow=false(size(FileList,1),1);
    
    for n=1:size(FileList,1)
        if FileList(n).isdir ~= 1
            if strcmpi(FileList(n).name(end-2:end),'tif')==1
                FileNames{n}=FileList(n).name;
                FileAllow(n)=true;
            end
        end
    end
    FileNames=FileNames(FileAllow);
    clear FileList FileAllow

NumExptPat=length(FileNames);

for countExptPat=1:NumExptPat

    Reference=imread([ExptPatFolder 'pattern_' num2str(countExptPat,'%5.5u') '.tif']);
    ReferenceImage= flipud(imresize(Reference,[H,W]));
% crop the data to remove the circle
    if crop_flag==1
        ReferenceImage = ReferenceImage(7:56,7:94);
    end

    % Loop through the simulation library
    for count=1:count2
        
        ImageToCompare=LibPattern{count};
        if crop_flag==1
            ImageToCompare = ImageToCompare(7:56,7:94);
        end
        % For each index as defined by the count, the correlation coefficient
        A{countExptPat}(count) = corr2(ReferenceImage, ImageToCompare);
        
    end


%%

    [CCvals{countExptPat},CCindex{countExptPat}] = sort(A{countExptPat},'descend');
    LibNum(countExptPat) = CCindex{countExptPat}(1);
    BestFit_CCval(countExptPat)=CCvals{countExptPat}(1);
    BestFit_phi1(countExptPat)=phi1(LibNum(countExptPat));
    BestFit_Phi(countExptPat)=Phi(LibNum(countExptPat));
    BestFit_phi2(countExptPat)=phi2(LibNum(countExptPat));
    testOri(countExptPat) = orientation('Euler',BestFit_phi1(countExptPat)*degree,BestFit_Phi(countExptPat)*degree,BestFit_phi2(countExptPat)*degree,CS,SS);
    BestFitOri(countExptPat)=testOri(countExptPat);
    MatchOut(countExptPat,:)=[countExptPat, LibNum(countExptPat), BestFit_CCval(countExptPat), BestFit_phi1(countExptPat), BestFit_Phi(countExptPat), BestFit_phi2(countExptPat)];


figure(101);
subplot(4,5,1);imagesc(imresize(ReferenceImage,[72,100]));axis image off xy;colormap('gray');
text(10,10,['pattern: ' num2str(countExptPat) '   CC: ' num2str(CCvals{countExptPat}(1)) ],'Color','red');
for i=1:19
subplot(4,5,i+1);imagesc(LibPattern{CCindex{countExptPat}(i)});axis image off xy;colormap('gray');
end
    saveas(101,[OutputFolder 'matches for Expt Pattern ' num2str(countExptPat) '.tif']);
    saveas(101,[OutputFolder 'matches for Expt Pattern ' num2str(countExptPat) '.png']);
    saveas(101,[OutputFolder 'matches for Expt Pattern ' num2str(countExptPat) '.fig']);

end


h = [Miller(1,0,0, CS)];
figure;
for countExptPat=1:NumExptPat
plotPDF(testOri(countExptPat),CCvals{countExptPat}(1),h,'projection','stereo','antipodal','MarkerSize',5, 'points', 1);colormap('jet');hold on;
end
    saveas(gcf,[OutputFolder '100_pole_figure.tif']);
    saveas(gcf,[OutputFolder '100_pole_figure.fig']);

h = [Miller(1,1,0, CS)];
figure;
for countExptPat=1:NumExptPat
plotPDF(testOri(countExptPat),CCvals{countExptPat}(1),h,'projection','stereo','antipodal','MarkerSize',5, 'points', 1);colormap('jet');hold on;
end
    saveas(gcf,[OutputFolder '110_pole_figure.tif']);
    saveas(gcf,[OutputFolder '110_pole_figure.fig']);

    elapsedTime = toc;
    hours=floor(elapsedTime/(60*60));
    mins=floor((elapsedTime-60*60*hours)/(60));
    secs=elapsedTime-60*60*hours - 60*mins;
    display(['pattern matching completed: hours=' num2str(hours) '  mins=' num2str(mins) '  secs=' num2str(secs)]);

MSAdata=load([InputFolder GrainFile '.mat']);
GrainMap=MSAdata.grainMap;    
for i=1:size(MSAdata.grainMap,1)
    for j=1:size(MSAdata.grainMap,2)
        OriNum=MSAdata.grainMap(i,j);
        CCvalMap(i,j)=CCvals{OriNum}(1);
    end
end
fname=[OutputFolder 'MatchMapData.mat'];
save(fname,'CCvalMap', 'BestFitOri', 'GrainMap' );   
    
    

    %%

oM = ipdfHSVOrientationMapping(CS);
oM.inversePoleFigureDirection = zvector;

figure;imagesc(MSAdata.grainMap);axis image;colormap('jet');colorbar;


tic
oM.inversePoleFigureDirection = zvector;
IPF_Z_color = oM.orientation2color(BestFitOri);
oM.inversePoleFigureDirection = yvector;
IPF_Y_color = oM.orientation2color(BestFitOri);
oM.inversePoleFigureDirection = xvector;
IPF_X_color = oM.orientation2color(BestFitOri);

Lcount=1;
fEBSD=fopen('tempEBSDmap.txt', 'w');
for i=1:size(MSAdata.grainMap,1)
    for j=1:size(MSAdata.grainMap,2)
        OriNum=MSAdata.grainMap(i,j);
        CCvalMap(i,j)=CCvals{OriNum}(1);
        LibNumMap(i,j)=CCindex{OriNum}(1);
        LibOriMap(i,j)=BestFitOri(OriNum);
        IPF_Z_map(i,j,1:3)=IPF_Z_color(OriNum,1:3);
        IPF_Y_map(i,j,1:3)=IPF_Y_color(OriNum,1:3);
        IPF_X_map(i,j,1:3)=IPF_X_color(OriNum,1:3);
%         fwrite(fEBSD,[i,j,OriNum,phi1(CCindex{OriNum}(1)),Phi(CCindex{OriNum}(1)),phi2(CCindex{OriNum}(1)),CCvalMap(i,j)]);
        EBSDtext(Lcount,:)=[i,j,OriNum,phi1(CCindex{OriNum}(1)),Phi(CCindex{OriNum}(1)),phi2(CCindex{OriNum}(1)),CCvalMap(i,j)];
        Lcount=Lcount+1;
    end
end
csvwrite('tempEBSDmap.csv',EBSDtext) ;
fclose(fEBSD);
toc

figure;imagesc(CCvalMap);axis image;colormap('jet');colorbar;


figure;image(IPF_Z_map);axis image off;title('IPF Z');
    saveas(gcf,[OutputFolder 'IPF Z map.png']);
    saveas(gcf,[OutputFolder 'IPF Z map.fig']);
figure;image(IPF_Y_map);axis image off;title('IPF Y');
    saveas(gcf,[OutputFolder 'IPF Y map.png']);
    saveas(gcf,[OutputFolder 'IPF Y map.fig']);
figure;image(IPF_X_map);axis image off;title('IPF X');
    saveas(gcf,[OutputFolder 'IPF X map.png']);
    saveas(gcf,[OutputFolder 'IPF X map.fig']);

    
figure;plot(oM);
    saveas(gcf,[OutputFolder 'IPF Colour Key.png']);
    saveas(gcf,[OutputFolder 'IPF Colour Key.fig']);
    
    
    
% for i=1:10
% imwrite(flipud(LibPattern{CCindex{2}(i)}), [OutputFolder 'Pattern2_fit' num2str(i) '.png'])
%     
% end

    
    
    
