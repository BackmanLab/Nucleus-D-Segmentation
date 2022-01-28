function [CellMarkerArray,AnalysisName,MaskfullName,MarkersName,MarkersValue]=NUPWS_D_CellLine126Markers_ExcelOutput(D_or_Sigma)
 
 %If "D_or_Sigma=1" then it will calculate the subsection results for D otherwise it does the analysis based on Sigma   
clc;
close all;  % Close all figures (except those of imtool.)
imtool close all;
clearvars -except D_or_Sigma  % Make sure the workspace panel is showing.
 
fontSize=22;
 
 
% Check if user has installed the Image Processing Toolbox.
hasIPT = license('test', 'image_toolbox');
if ~hasIPT
  % the Image Processing Toolbox is not installed.
  message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
  reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
  if strcmpi(reply, 'No')
    % User said No, so exit.
    return;
  end
end
 
%[baseFileName,folder] = uigetfile
%fullFileName = fullfile(folder, baseFileName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Have user browse for a file, from a specified "starting folder."
% For convenience in browsing, set a starting folder from which to browse.
startingFolder = 'C:\Program Files\MATLAB';
if ~exist(startingFolder, 'dir')
  % If that folder doesn't exist, just start in the current folder.
  startingFolder = pwd;
end
% Get the name of the file that the user wants to use.
defaultFileName = fullfile(startingFolder, '*.*');
[baseFileName, folder] = uigetfile(defaultFileName, 'Select a file');
if baseFileName == 0
  % User clicked the Cancel button.
  return;
end
fullFileName = fullfile(folder, baseFileName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=1;
for v=1:15
    v
clearvars -except D_or_Sigma CellMarkerArray MaskRow DatasetRow MaskfullName AnalysisName v fullFileName k folder baseFileName kk
 
idx1 = strfind(fullFileName,"Cell");
idx2 = strfind(fullFileName,"PWS");
fullFileName(idx1+4:idx2-2)=[];
CellNumber=num2str(v);
fullFileName=insertAfter(fullFileName,"Cell",CellNumber);
 
 
idx1 = strfind(folder,"Cell");
idx2 = strfind(folder,"PWS");
folder(idx1+4:idx2-2)=[];
folder=insertAfter(folder,"Cell",CellNumber);
 
 
 
%a = h5read('C:\Users\ali_n\Documents\Northwestern PC\Projects\PWS Project\Projects\SOP\95% Ethanol Fixation _Time Variant Analysis\Hella Cell_Practice\liveorig\Cell E\Cell52\PWS\analyses\analysisResults_fixed.h5', '/rms')
NW_PWSImage = h5read(fullFileName, '/rms');
noise=0.0372;
true_rms = sqrt(abs(NW_PWSImage.^2 - noise.^2));
NW_PWSImage=true_rms;

f=folder;
f1=folder(1:end-13);
d=dir(f1);
d_folder=d.folder;
dd=d(end);
kk=1;
 
for ii=1:3 %Find ROI 01,   ROI 02 and ROI 03.
    dd=d(end+1-ii);
    hh=dd.name(end-1:end);
if hh=='h5'
    maskname{ii}=dd.name;
    maskfolder{ii}=dd.folder;
    maskfullFileName{ii}=fullfile(dd.folder, dd.name);
 
    info=h5info(maskfullFileName{ii});
 
    %info.Groups.Name
    GroupNames=struct2cell(info.Groups);
    [m n]=size(GroupNames);
    
for j=1:n
%GroupName=struct2cell(info.Groups.Name)
GroupName=char(GroupNames(1,j));
DatasetName=char({'/mask/'});
as=([GroupName DatasetName]);
msk = h5read(maskfullFileName{ii}, as);
%m4 = h5read(maskfullFileName{1,1}, '/4/mask/');
Mask(ii,j,:,:)=msk;
 
MaskRow(k,:,:)=msk;
DatasetRow{k}=GroupName;
MaskfullName{k}=fullfile(dd.folder, dd.name);
AnalysisName{k}=f;
 
Mask_Row(v,kk,:,:)=msk;
Dataset_Row{v,kk}=GroupName;
Mask_fullName{v,kk}=fullfile(dd.folder, dd.name);
Analysis_Name{v,kk}=f;
 
 
k=k+1;
kk=kk+1;
NucMask=msk;
 
binaryImage = imbinarize(10*NW_PWSImage);
%imshow(binaryImage)
BWfilter = bwpropfilt(binaryImage,'perimeter',5);
%imshow(BWfilter)
 
 
windowSize = 16; % Whatever you want.
kernel = ones(windowSize, windowSize) / windowSize ^ 2;
blurredBWfilter = imfilter(BWfilter, kernel, 'symmetric');
CellMask= bwpropfilt(blurredBWfilter,'perimeter',5);
CellImage=NW_PWSImage;

CellMask(1:45,980:1024)=0;
CellMask(1:45,1:45)=0;
CellMask(980:1024,1:45)=0;
CellMask(980:1024,980:1024)=0;
CellImage(CellMask==0)=0;

if D_or_Sigma==1
%Transfer Sigma Map to Dmap
phi=0.35;Nf=1e6;thickness=2;Sigma_min=0.02;Sigma_max=0.6;
%Phi is estimated based on A549 cells (between 0.3-0.4).Thickness was
%measured for HCT116 to be 2 micron by 3D optical profilometer.
[polyVals]=NU_SigmaToD_polyApprox(phi,Nf,thickness,Sigma_min,Sigma_max);
% Polynomial function that estimates D from Sigma values
 D_map= polyval(polyVals, CellImage);
 CellImage=D_map;
else
  CellImage=CellImage;
end
%Random change


[MakerArray]=PWSLocalD126Markers(CellImage,NucMask);
CellMarkerArray(k,:)=MakerArray;
Cell_Marker_Array(v,kk,:)=MakerArray;
 
end
 
 
end
end
end
MaskfullName=MaskfullName';
DatasetRow=DatasetRow';
AnalysisName=AnalysisName';
 
[nn mm]=size(CellMarkerArray);
 
 

f1=folder;f1(end-19:end)=[];
f2='_D_126markersMatlabData_Aug25_2021';
f3='.xlsx';
f=[f1 f2 f3];

AnalysisName = erase(AnalysisName,f1);

xlsfilename=f; 
MarkersName={'Marker1';'Marker2';'Marker3';'Marker4';'Marker5';'Marker6';'Marker7';'Marker8';'Marker9';'Marker10';'Marker11';'Marker12';'Marker13';'Marker14';'Marker15';'Marker16';'Marker17';'Marker18';'Marker19';'Marker20';'Marker21';'Marker22';'Marker23';'Marker24';'Marker25';'Marker26';'Marker27';'Marker28';'Marker29';'Marker30';'Marker31';'Marker32';'Marker33';'Marker34';'Marker35';'Marker36';'Marker37';'Marker38';'Marker39';'Marker40';'Marker41';'Marker42';'Marker43';'Marker44';'Marker45';'Marker46';'Marker47';'Marker48';'Marker49';'Marker50';'Marker51';'Marker52';'Marker53';'Marker54';'Marker55';'Marker56';'Marker57';'Marker58';'Marker59';'Marker60';'Marker61';'Marker62';'Marker63';'Marker64';'Marker65';'Marker66';'Marker67';'Marker68';'Marker69';'Marker70';'Marker71';'Marker72';'Marker73';'Marker74';'Marker75';'Marker76';'Marker77';'Marker78';'Marker79';'Marker80';'Marker81';'Marker82';'Marker83';'Marker84';'Marker85';'Marker86';'Marker87';'Marker88';'Marker89';'Marker90';'Marker91';'Marker92';'Marker93';'Marker94';'Marker95';'Marker96';'Marker97';'Marker98';'Marker99';'Marker100';'Marker101';'Marker102';'Marker103';'Marker104';'Marker105';'Marker106';'Marker107';'Marker108';'Marker109';'Marker110';'Marker111';'Marker112';'Marker113';'Marker114';'Marker115';'Marker116';'Marker117';'Marker118';'Marker119';'Marker120';'Marker121';'Marker122';'Marker123';'Marker124';'Marker125';'Marker126'}';
MarkersName=string(MarkersName);

% Subsection Type 1: Average D,AreaRatio,ConnectedArea1,Std,ConnectedArea2,ConnectedArea3 Variables in each subsection
%Markesr# 1 to 48
AveRms_Subsectopn1={'AveRms (<µCell-2S)';'AveRms (<µCell-S)';'AveRms (<µCell)';'AveRms (>µCell)';'AveRms (>µCell+S)';'AveRms (>µCell+2S)';'AveRms (>µCell+3S)';'AveRms (>µCell+4S)'}
AreaRatio_Subsectopn1={'AreaRatio (<µCell-2S)';'AreaRatio (<µCell-S)';'AreaRatio (<µCell)';'AreaRatio (>µCell)';'AreaRatio (>µCell+S)';'AreaRatio (>µCell+2S)';'AreaRatio (>µCell+3S)';'AreaRatio (>µCell+4S)'}
ConnectedArea1_Subsectopn1={'ConnectedArea1 (<µCell-2S)';'ConnectedArea1 (<µCell-S)';'ConnectedArea1 (<µCell)';'ConnectedArea1 (>µCell)';'ConnectedArea1 (>µCell+S)';'ConnectedArea1 (>µCell+2S)';'ConnectedArea1 (>µCell+3S)';'ConnectedArea1 (>µCell+4S)'}
StdRms_Subsectopn1={'StdRms (<µCell-2S)';'StdRms (<µCell-S)';'StdRms (<µCell)';'StdRms (>µCell)';'StdRms (>µCell+S)';'StdRms (>µCell+2S)';'StdRms (>µCell+3S)';'StdRms (>µCell+4S)';}
ConnectedArea2_Subsectopn1={'ConnectedArea2 (<µCell-2S)';'ConnectedArea2 (<µCell-S)';'ConnectedArea2 (<µCell)';'ConnectedArea2 (>µCell)';'ConnectedArea2 (>µCell+S)';'ConnectedArea2 (>µCell+2S)';'ConnectedArea2 (>µCell+3S)';'ConnectedArea2 (>µCell+4S)'}
ConnectedArea3_Subsectopn1={'ConnectedArea3 (<µCell-2S)';'ConnectedArea3 (<µCell-S)';'ConnectedArea3 (<µCell)';'ConnectedArea3 (>µCell)';'ConnectedArea3 (>µCell+S)';'ConnectedArea3 (>µCell+2S)';'ConnectedArea3 (>µCell+3S)';'ConnectedArea3 (>µCell+4S)'}
MarkersValue_Sucsection1=[AveRms_Subsectopn1;AreaRatio_Subsectopn1;ConnectedArea1_Subsectopn1;StdRms_Subsectopn1;ConnectedArea2_Subsectopn1;ConnectedArea3_Subsectopn1]

%Marker#49
TraditionalRms={'Traditional Ave Rms'};

% Subsection Type 1 (Lables): Average D,AreaRatio,ConnectedArea1,Std,ConnectedArea2,ConnectedArea3 Variables in each subsection
%%Marker# 49 to 91
AveRms_Subsectopn2={'µCell-2S<AveRms<µCell-S';'µCell-S<AveRms<µCell';'µCell<AveRms<µCell+S';'µCell+S<AveRms<µCell+2S';'µCell+2S<AveRms<µCell+3S';'µCell+3S<AveRms<µCell+4S';'AveRms (>µCell+4S)'}
AreaRatio_Subsectopn2={'µCell-2S<AreaRatio<µCell-S';'µCell-S<AreaRatio<µCell';'µCell<AreaRatio<µCell+S';'µCell+S<AreaRatio<µCell+2S';'µCell+2S<AreaRatio<µCell+3S';'µCell+3S<AreaRatio<µCell+4S';'AreaRatio (>µCell+4S)'}
ConnectedArea1_Subsectopn2={'µCell-2S<ConnectedArea1<µCell-S';'µCell-S<ConnectedArea1<µCell';'µCell<ConnectedArea1<µCell+S';'µCell+S<ConnectedArea1<µCell+2S';'µCell+2S<ConnectedArea1<µCell+3S';'µCell+3S<ConnectedArea1<µCell+4S';'ConnectedArea1 (>µCell+4S)'}
StdRms_Subsectopn2={'µCell-2S<StdRms<µCell-S';'µCell-S<StdRms<µCell';'µCell<StdRms<µCell+S';'µCell+S<StdRms<µCell+2S';'µCell+2S<StdRms<µCell+3S';'µCell+3S<StdRms<µCell+4S';'StdRms (>µCell+4S)'}
ConnectedArea2_Subsectopn2={'µCell-2S<ConnectedArea2<µCell-S';'µCell-S<ConnectedArea2<µCell';'µCell<ConnectedArea2<µCell+S';'µCell+S<ConnectedArea2<µCell+2S';'µCell+2S<ConnectedArea2<µCell+3S';'µCell+3S<ConnectedArea2<µCell+4S';'ConnectedArea2 (>µCell+4S)'}
ConnectedArea3_Subsectopn2={'µCell-2S<ConnectedArea3<µCell-S';'µCell-S<ConnectedArea3<µCell';'µCell<ConnectedArea3<µCell+S';'µCell+S<ConnectedArea3<µCell+2S';'µCell+2S<ConnectedArea3<µCell+3S';'µCell+3S<ConnectedArea3<µCell+4S';'ConnectedArea3 (>µCell+4S)'}
MarkersValue_Sucsection2=[AveRms_Subsectopn2;AreaRatio_Subsectopn2;ConnectedArea1_Subsectopn2;StdRms_Subsectopn2;ConnectedArea2_Subsectopn2;ConnectedArea3_Subsectopn2];

%Differantial Markers (Lables) 92 to 111
Differential_DMarkers={'PD(D>µCell)-PD(D<µCell)';'PD(D<µCell)-PD(D<µCell-S)';'PD(D<µCell-S)-PD(D<µCell-2S)';'PD(D>µCell+S-PD(D>µCell)';'PD(D>µCell+2S-PD(D>µCell+S)';'PD(D>µCell+3S)-PD(D>µCell+2S)';'PD(D>µCell+4S)-PD(D>µCell+3S)'};
Differential_AreaMarkers={'Area(D<µCell-S)-Area(D<µCell-2S)';'Area(D<µCell)-PD(D<µCell-S)';'Area(D>µCell)-PD(D<µCell)';'Area(D>µCell+S-Area(D>µCell)';'Area(D>µCell+2S-Area(D>µCell+S)';'Area(D>µCell+3S)-Area(D>µCell+2S)';'Area(D>µCell+4S)-Area(D>µCell+3S)'};
Differential_StdMarkers={'LowD D-DStd(Markers 92,93,94)';'HighD D-DStd(Markers 95,96,97,98)';'All D-DStd(Markers 92to98)';'LowD D-AreaStd(Markers 99,100,101)';'HighD D-AreaStd(Markers 102,103,104,105)';'All D-Area-Std(Markers 92 to 97)'  }
DifferentialMarkers=[Differential_DMarkers;Differential_AreaMarkers;Differential_StdMarkers];

%Boundry Markers (Lables) Markers 112 to 123
Boundry_Markers={'Boundry_Edge_Ave';'Boundry_NearEdge_Ave';'Boundry_Center_Ave';'Boundry_Edge_Std';'Boundry_NearEdge_Std';'Center_Std';'Boundry_Edge_entropy';'Boundry_NearEdge_entropy';'Center_entropy';'Edge/NearEdge_Ave';'Edge/Center_Ave';'NearEdge/Center_Ave';'Edge*NearEdge_Ave';'Edge*Center';'NearEdge*Center_Ave'}
%                 Boundry_Edge_Av Boundry_NearEdge_Ave	Center_Ave	            Boundry_Edge_Std	Boundry_NearEdge_Std	Center_Std	Boundry_Edge_entropy	Boundry_NearEdge_entropy	Center_entropy	Edge/NearEdge_Ave Edge/Center_Ave	NearEdge/Center_Ave	   Edge*NearEdge_Ave	Edge*Center	NearEdge*Center_Ave
MarkersValue=[MarkersValue_Sucsection1;TraditionalRms;MarkersValue_Sucsection2;DifferentialMarkers;Boundry_Markers];

MarkersValue=string(MarkersValue);

MarkersValue=MarkersValue';

ii=1;
firstRow =(ii-1)*60+5;
lastRow = (ii-1)*60+5+nn;
firstCol = 'D';
lastCol = 'DY';

cellRange = [firstCol,num2str(firstRow),':',lastCol,num2str(lastRow)];
Header1Range=[firstCol,num2str(firstRow-1),':',lastCol,num2str(lastRow-1)];
Header2Range=['A',num2str(firstRow)+1,':','D',num2str(lastRow+1)];
Header3Range=['B',num2str(firstRow)+1,':','D',num2str(lastRow+1)];
Header4Range=['C',num2str(firstRow)+1,':','D',num2str(lastRow+1)];
Header5Range=[firstCol,num2str(firstRow),':',lastCol,num2str(lastRow)];
Header6Range=[firstCol,num2str(firstRow-1),':',lastCol,num2str(lastRow-1)];

 
writematrix(CellMarkerArray,xlsfilename,'Sheet',1,'Range',cellRange)
writecell(DatasetRow,xlsfilename,'Sheet',1,'Range',Header2Range)
writecell(AnalysisName,xlsfilename,'Sheet',1,'Range',Header3Range)
writecell(MaskfullName,xlsfilename,'Sheet',1,'Range',Header4Range)
writematrix(MarkersName,xlsfilename,'Sheet',1,'Range',Header5Range)
writematrix(MarkersValue,xlsfilename,'Sheet',1,'Range',Header6Range)



%Ali Daneshkhah, PhD 
%Postdoctoral Fellow
%Backman Photonics Laboratory 
%Biomedical Engineering Department
%McCormick School of Engineering & Applied Sciences
%Northwestern University