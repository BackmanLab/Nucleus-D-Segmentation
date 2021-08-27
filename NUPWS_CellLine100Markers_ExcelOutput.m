function [CellMarkerArray,AnalysisName,MaskfullName,MarkersName,MarkersValue]=NUPWS_CellLine100Markers_ExcelOutput()
 
    
clc;
close all;  % Close all figures (except those of imtool.)
imtool close all;
clear all
workspace;  % Make sure the workspace panel is showing.
 
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
 
[baseFileName,folder] = uigetfile
fullFileName = fullfile(folder, baseFileName);
k=1;
for v=1:4
    %v
clearvars -except CellMarkerArray MaskRow DatasetRow MaskfullName AnalysisName v fullFileName k folder baseFileName kk
 
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

AveRMS=mean(mean(NW_PWSImage));
 
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
CC(CellMask==0)=0;

[MakerArray]=PWSLocal100Markers(CellImage,NucMask);
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
f2='_Rms100markersMatlabData_June29_2021';
f3='.xlsx';
f=[f1 f2 f3];

AnalysisName = erase(AnalysisName,f1);

xlsfilename=f; 
MarkersName={'Marker1';'Marker2';'Marker3';'Marker4';'Marker5';'Marker6';'Marker7';'Marker8';'Marker9';'Marker10';'Marker11';'Marker12';'Marker13';'Marker14';'Marker15';'Marker16';'Marker17';'Marker18';'Marker19';'Marker20';'Marker21';'Marker22';'Marker23';'Marker24';'Marker25';'Marker26';'Marker27';'Marker28';'Marker29';'Marker30';'Marker31';'Marker32';'Marker33';'Marker34';'Marker35';'Marker36';'Marker37';'Marker38';'Marker39';'Marker40';'Marker41';'Marker42';'Marker43';'Marker44';'Marker45';'Marker46';'Marker47';'Marker48';'Marker49';'Marker50';'Marker51';'Marker52';'Marker53';'Marker54';'Marker55';'Marker56';'Marker57';'Marker58';'Marker59';'Marker60';'Marker61';'Marker62';'Marker63';'Marker64';'Marker65';'Marker66';'Marker67';'Marker68';'Marker69';'Marker70';'Marker71';'Marker72';'Marker73';'Marker74';'Marker75';'Marker76';'Marker77';'Marker78';'Marker79';'Marker80';'Marker81';'Marker82';'Marker83';'Marker84';'Marker85';'Marker86';'Marker87';'Marker88';'Marker89';'Marker90';'Marker91';'Marker92';'Marker93';'Marker94';'Marker95';'Marker96';'Marker97';'Marker98';'Marker99';'Marker100'}';
MarkersName=string(MarkersName);

MarkersValue={'AveRms (>�Cell)';'Area Ratio (>�Cell)';'Rms Std (>�Cell)';'AveRms (�<x<�+S)';'Area Ratio (�<x<�+? )';'RmsStd (�<x<�+?)';'AveRms (>�Cell+S)';'Area Ratio (>�Cell+?)';'Rms Std (>�Cell+?)';'AveRms (�+?<x<�+2?)';'Area Ratio (�+?<x<�+2? )';'RmsStd (�+?<x<�+2?)';'AveRms (>�Cell+2?)';'Area Ratio (>�Cell+2?)';'Rms Std (>�Cell+2?)';'AveRms (�+2?<x<�+3?)';'Area Ratio (�+2?<x<�+3? )';'RmsStd (�+2?<x<�+3?)';'AveRms (>�Cell+3?)';'Area Ratio (>�Cell+3?)';'Rms Std (>�Cell+3?)';'AveRms (�+3?<x<�+4?)';'Area Ratio (�+3?<x<�+4? )';'RmsStd (�+3?<x<�+4?)';'AveRms (>�Cell+4?)';'Area Ratio (>�Cell+4?)';'Rms Std (>�Cell+4?)';'AveRms (>�Nuc)';'Area Ratio (>�Nucl)';'Rms Std (>�Nuc)';'AveRms (>�Cell+?)-AveRms (>�Cell)';'AveArea (>�Cell+?)-AveRms (>�Cell)';'AveRms (>�Cell+2?)-AveRms (>�Cell+?)';'AveArea (>�Cell+2?)-AveRms (>�Cell+?)';'AveRms (>�Cell+3?)-AveRms (>�Cell+2?)';'AveArea (>�Cell+3?)-AveRms (>�Cell+2?)';'AveRms (>�Cell+4?)-AveRms (>�Cell+3?)';'AveArea (>�Cell+4?)-AveArea (>�Cell+3?)';'Std (Markers_31,33,35,37)/4';'Std (Markers_32,34,36,58)/4';'Ave( Markers_31,33,35,37)';'Ave( Markers_32,34,36,37)';'Ave_Rms';	'AveRms (<�Cell)';'Area Ratio (<�Cell)';'Rms Std (<�Cell)';'AveRms (�-?<x<� )';'Area Ratio (�-?<x<� )';'RmsStd (�-?<x<� )';'AveRms (x<�-? )';'Area Ratio (x<�-? )';'RmsStd (x<�-? )';'AveRms (�-2?<x<�-? )';'Area Ratio (�-2?<x<�-? )';'RmsStd (�-2?<x<� )';'ConnectedArea1 (>�Cell)';'ConnectedArea2 (>�Cell)';'ConnectedArea3 (>�Cell)';'ConnectedArea1 (�<x<�+S)';'ConnectedArea2 (�<x<�+S)';'ConnectedArea3 (�<x<�+S)';'ConnectedArea1 (>�Cell+S)';'ConnectedArea2 (>�Cell+S)';'ConnectedArea3 (>�Cell+S)';'ConnectedArea1 (�+S<x<�+2S)';'ConnectedArea2 (�+S<x<�+2S)';'ConnectedArea3 (�+S<x<�+2S)';'ConnectedArea1 (>�Cell+2S)';'ConnectedArea2 (>�Cell+2S)';'ConnectedArea3 (>�Cell+2S)';'ConnectedArea1 (�+2S<x<�+3S)';'ConnectedArea2 (�+2S<x<�+3S)';'ConnectedArea3 (�+2S<x<�+3S)';'ConnectedArea1 (>�Cell+3S)';'ConnectedArea2 (>�Cell+3S)';'ConnectedArea3 (>�Cell+3S)';'ConnectedArea1 (�+3S<x<�+4S)';'ConnectedArea2 (�+3S<x<�+4S)';'ConnectedArea3 (�+3S<x<�+4S)';'ConnectedArea1 (>�Cell+4S)';'ConnectedArea2 (>�Cell+4S)';'ConnectedArea3 (>�Cell+4S)';'ConnectedArea1 (<�Cell)';'ConnectedArea2 (<�Cell)';'ConnectedArea3 (<�Cell)';'ConnectedArea1 (�-?<x<� )';'ConnectedArea2 (�-?<x<� )';'ConnectedArea3 (�-?<x<� )';'ConnectedArea1 (x<�-? )';'ConnectedArea2 (x<�-? )';'ConnectedArea3 (x<�-? )';'ConnectedArea1 (�-2?<x<� )';'ConnectedArea2 (�-2?<x<� )';'ConnectedArea3 (�-2?<x<� )';'Marker_Boundry1';'Marker_Boundry2';'Marker_Center';'Marker_Boundry1/Marker_Boundry2';'Marker_Boundry1/Marker_Center';'Marker_Boundry2/Marker_Center'}';
MarkersValue=string(MarkersValue);


ii=1;
firstRow =(ii-1)*60+5;
lastRow = (ii-1)*60+5+nn;
firstCol = 'D';
lastCol = 'CY';

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