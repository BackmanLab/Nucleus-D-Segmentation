clear all
close all
clc
[baseFileName,folder] = uigetfile


clearvars -except baseFileName folder
for ii=46:46
folder='/Volumes/My Passport/PWS_Data_Liquid/BMC2018/#546_2/Cell1/';
%folder='/Volumes/My Passport/PWS_Data_Liquid/NMH-DM2017/DM_Control/#3_2/Cell1/'
%baseFileName='450_648_f6_30_gp2_adc_9981_Ld.mat'
total_cells=39;

    
str0 = append('#',num2str(546));
newstr0=append('#',num2str(555+ii));
folder = strrep(folder,str0,newstr0);
    
    
for i=1:total_cells
clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
close all

if i>1
str = append('Cell',num2str(i-1));
newstr=append('Cell',num2str(i));
folder = strrep(folder,str,newstr);
end


[NCOutput NCOutputCell NCOutputNuc maskNuc Marker_address]=NCImageOutputRev1b(baseFileName,folder);
ImageCell(:,:)=NCOutputCell(4,:,:);
%Mask(:,:)=NCOutput(18,:,:);
Mask=maskNuc;
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(2)

All_Cell_Image_Array4(i,:,:)=MakerArray;
close all
clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(12,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array12(i,:,:)=MakerArray;

%%%%%PWSLocal100Markers

clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(1,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array1(i,:,:)=MakerArray;

%%%%%

clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(1,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array1(i,:,:)=MakerArray;

%%%%%


clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(2,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array2(i,:,:)=MakerArray;


%%%%%

clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(3,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array3(i,:,:)=MakerArray;

%%%%%

%%%%%


clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(5,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array5(i,:,:)=MakerArray;


clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(6,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array6(i,:,:)=MakerArray;

clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(7,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array7(i,:,:)=MakerArray;

clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(8,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array8(i,:,:)=MakerArray;

clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(9,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array9(i,:,:)=MakerArray;

clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(10,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array10(i,:,:)=MakerArray;


clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(11,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array11(i,:,:)=MakerArray;


clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(13,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array13(i,:,:)=MakerArray;

clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(14,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array14(i,:,:)=MakerArray;

clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(15,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array15(i,:,:)=MakerArray;


clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16 S2 ii Marker_address All_Cell_Image_Array1 All_Cell_Image_Array2 All_Cell_Image_Array3 All_Cell_Image_Array16
ImageCell(:,:)=NCOutputCell(16,:,:);
[MakerArray]=PWSLocal100Markers(ImageCell,Mask);
%pause(1)
All_Cell_Image_Array16(i,:,:)=MakerArray;

%clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder total_cells All_Cell_Image_Array5 All_Cell_Image_Array6 All_Cell_Image_Array7 All_Cell_Image_Array8 All_Cell_Image_Array9 All_Cell_Image_Array10 All_Cell_Image_Array11 All_Cell_Image_Array13 All_Cell_Image_Array14 All_Cell_Image_Array15 All_Cell_Image_Array16
%ImageCell(:,:)=NCOutputCell(16,:,:);
%[MakerArray]=PWSLocal43Markers(ImageCell,Mask);
%pause(1)
%All_Cell_Image_Array16(i,:,:)=MakerArray;

%Pause(1)
close all
%clearvars -except All_Cell_Image_Array4 All_Cell_Image_Array12 All_Cell_Image_Array14 i NCOutput NCOutputCell NCOutputNuc Mask baseFileName folder
%ImageCell(:,:)=NCOutputCell(14,:,:);
%[MakerArray]=PWSLocalMarkers(ImageCell,Mask);
%pause(1)
%All_Cell_Image_Array14(i,:,:)=MakerArray;


i
end

a1(:,:)=All_Cell_Image_Array1(1:total_cells,1,:);
a2(:,:)=All_Cell_Image_Array2(1:total_cells,1,:);
a3(:,:)=All_Cell_Image_Array3(1:total_cells,1,:);

a4(:,:)=All_Cell_Image_Array4(1:total_cells,1,:);
a12(:,:)=All_Cell_Image_Array12(1:total_cells,1,:);


a5(:,:)=All_Cell_Image_Array5(1:total_cells,1,:);
a6(:,:)=All_Cell_Image_Array6(1:total_cells,1,:);
a7(:,:)=All_Cell_Image_Array7(1:total_cells,1,:);
a8(:,:)=All_Cell_Image_Array8(1:total_cells,1,:);
a9(:,:)=All_Cell_Image_Array9(1:total_cells,1,:);
a10(:,:)=All_Cell_Image_Array10(1:total_cells,1,:);
a11(:,:)=All_Cell_Image_Array11(1:total_cells,1,:);
a13(:,:)=All_Cell_Image_Array13(1:total_cells,1,:);
a14(:,:)=All_Cell_Image_Array14(1:total_cells,1,:);
a15(:,:)=All_Cell_Image_Array15(1:total_cells,1,:);
a16(:,:)=All_Cell_Image_Array16(1:total_cells,1,:);


%a14(:,:)=All_Cell_Image_Array14(1:17,1,:);
%Questions can be adressed to alid@northwestern.edu








%%%%----------------

%ii=5;
MarkersName={'Marker1';'Marker2';'Marker3';'Marker4';'Marker5';'Marker6';'Marker7';'Marker8';'Marker9';'Marker10';'Marker11';'Marker12';'Marker13';'Marker14';'Marker15';'Marker16';'Marker17';'Marker18';'Marker19';'Marker20';'Marker21';'Marker22';'Marker23';'Marker24';'Marker25';'Marker26';'Marker27';'Marker28';'Marker29';'Marker30';'Marker31';'Marker32';'Marker33';'Marker34';'Marker35';'Marker36';'Marker37';'Marker38';'Marker39';'Marker40';'Marker41';'Marker42';'Marker43';'Marker44';'Marker45';'Marker46';'Marker47';'Marker48';'Marker49';'Marker50';'Marker51';'Marker52';'Marker53';'Marker54';'Marker55';'Marker56';'Marker57';'Marker58';'Marker59';'Marker60';'Marker61';'Marker62';'Marker63';'Marker64';'Marker65';'Marker66';'Marker67';'Marker68';'Marker69';'Marker70';'Marker71';'Marker72';'Marker73';'Marker74';'Marker75';'Marker76';'Marker77';'Marker78';'Marker79';'Marker80';'Marker81';'Marker82';'Marker83';'Marker84';'Marker85';'Marker86';'Marker87';'Marker88';'Marker89';'Marker90';'Marker91';'Marker92';'Marker93';'Marker94';'Marker95';'Marker96';'Marker97';'Marker98';'Marker99';'Marker100'}';
MarkersName=string(MarkersName);



f1=folder;f1(end-6:end)=[];
f2='_NCRms100markersMatlabData_July14_2021';
f3='.xlsx';
f=[f1 f2 f3];

%%%%%%%----------- Make a Master Excel file that includes all the Samples
%%%%%%%data
g='/Volumes/My Passport/PWS_Data_Liquid/BMC2018/#NC100Markers_July2021_.xlsx'
xlsfilename=g;

 firstRow =(ii-1)*60+3;
lastRow = (ii-1)*60+58;
firstCol = 'D';
lastCol = 'CY';
cellRange = [firstCol,num2str(firstRow),':',lastCol,num2str(lastRow)];
Header1Range=[firstCol,num2str(firstRow-1),':',lastCol,num2str(lastRow-1)];
Header2Range=[firstCol,num2str(firstRow-2)];
Header3Range=['E',num2str(firstRow-2)];
Header4Range=['F',num2str(firstRow-2)];


writematrix(a1,xlsfilename,'Sheet',1,'Range',cellRange)

writematrix(MarkersName,xlsfilename,'Sheet',1,'Range',Header1Range)
writematrix('A1_LD_gp0',xlsfilename,'Sheet',1,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',1,'Range',Header3Range)
writematrix(char(Marker_address(1)),xlsfilename,'Sheet',1,'Range',Header4Range)


writematrix(a2,xlsfilename,'Sheet',2,'Range',cellRange)

writematrix(MarkersName,xlsfilename,'Sheet',2,'Range',Header1Range)
writematrix('A2_Rsquared_gp0',xlsfilename,'Sheet',2,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',2,'Range',Header3Range)
writematrix(char(Marker_address(2)),xlsfilename,'Sheet',2,'Range',Header4Range)



writematrix(a3,xlsfilename,'Sheet',3,'Range',cellRange)

writematrix(MarkersName,xlsfilename,'Sheet',3,'Range',Header1Range)
writematrix('A3_Reflectance_gp0',xlsfilename,'Sheet',3,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',3,'Range',Header3Range)
writematrix(char(Marker_address(3)),xlsfilename,'Sheet',3,'Range',Header4Range)




writematrix(a4,xlsfilename,'Sheet',4,'Range',cellRange)

writematrix(MarkersName,xlsfilename,'Sheet',4,'Range',Header1Range)
writematrix('A4_RMS_gp0',xlsfilename,'Sheet',4,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',4,'Range',Header3Range)
writematrix(char(Marker_address(4)),xlsfilename,'Sheet',4,'Range',Header4Range)



writematrix(a5,xlsfilename,'Sheet',5,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',5,'Range',Header1Range)
writematrix('A5_RmsNorm_gp0',xlsfilename,'Sheet',5,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',5,'Range',Header3Range)
writematrix(char(Marker_address(5)),xlsfilename,'Sheet',5,'Range',Header4Range)


writematrix(a6,xlsfilename,'Sheet',6,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',6,'Range',Header1Range)
writematrix('A6_RmsPoly_gp0',xlsfilename,'Sheet',6,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',6,'Range',Header3Range)
writematrix(char(Marker_address(6)),xlsfilename,'Sheet',6,'Range',Header4Range)


writematrix(a7,xlsfilename,'Sheet',7,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',7,'Range',Header1Range)
writematrix('A7_RmsSquared_gp0',xlsfilename,'Sheet',7,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',7,'Range',Header3Range)
writematrix(char(Marker_address(7)),xlsfilename,'Sheet',7,'Range',Header4Range)


writematrix(a8,xlsfilename,'Sheet',8,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',8,'Range',Header1Range)
writematrix('A8_Slope_gp0',xlsfilename,'Sheet',8,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',8,'Range',Header3Range)
writematrix(char(Marker_address(8)),xlsfilename,'Sheet',8,'Range',Header4Range)


writematrix(a9,xlsfilename,'Sheet',9,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',9,'Range',Header1Range)
writematrix('A9_LD_gp2',xlsfilename,'Sheet',9,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',9,'Range',Header3Range)
writematrix(char(Marker_address(9)),xlsfilename,'Sheet',9,'Range',Header4Range)


writematrix(a10,xlsfilename,'Sheet',10,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',10,'Range',Header1Range)
writematrix('A10_RSquared_gp2',xlsfilename,'Sheet',10,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',10,'Range',Header3Range)
writematrix(char(Marker_address(10)),xlsfilename,'Sheet',10,'Range',Header4Range)


writematrix(a11,xlsfilename,'Sheet',11,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',11,'Range',Header1Range)
writematrix('A11_Reflectance_gp2',xlsfilename,'Sheet',11,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',11,'Range',Header3Range)
writematrix(char(Marker_address(11)),xlsfilename,'Sheet',11,'Range',Header4Range)


writematrix(a12,xlsfilename,'Sheet',12,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',12,'Range',Header1Range)
writematrix('A12_RMS_gp2',xlsfilename,'Sheet',12,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',12,'Range',Header3Range)
writematrix(char(Marker_address(12)),xlsfilename,'Sheet',12,'Range',Header4Range)


writematrix(a13,xlsfilename,'Sheet',13,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',13,'Range',Header1Range)
writematrix('A13_RmsNorm_gp2',xlsfilename,'Sheet',13,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',13,'Range',Header3Range)
writematrix(char(Marker_address(13)),xlsfilename,'Sheet',13,'Range',Header4Range)


writematrix(a14,xlsfilename,'Sheet',14,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',14,'Range',Header1Range)
writematrix('A14_RmsPoly_gp2',xlsfilename,'Sheet',14,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',14,'Range',Header3Range)
writematrix(char(Marker_address(14)),xlsfilename,'Sheet',14,'Range',Header4Range)



writematrix(a15,xlsfilename,'Sheet',15,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',15,'Range',Header1Range)
writematrix('A15_RmsSquared_gp2',xlsfilename,'Sheet',15,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',15,'Range',Header3Range)
writematrix(char(Marker_address(15)),xlsfilename,'Sheet',15,'Range',Header4Range)


writematrix(a16,xlsfilename,'Sheet',16,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',16,'Range',Header1Range)
writematrix('A16_Slope_gp2',xlsfilename,'Sheet',16,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',16,'Range',Header3Range)
writematrix(char(Marker_address(16)),xlsfilename,'Sheet',16,'Range',Header4Range)


%%%%%%%%%% 




%%%%%%%----------- Make a new Excel file for each Samples data
xlsfilename=f;
i=1;

 firstRow =(i-1)*50+3;
lastRow = (i-1)*50+55;
firstCol = 'D';
lastCol = 'CY';
cellRange = [firstCol,num2str(firstRow),':',lastCol,num2str(lastRow)];
Header1Range=[firstCol,num2str(firstRow-1),':',lastCol,num2str(lastRow-1)];
Header2Range=[firstCol,num2str(firstRow-2)];
Header3Range=['E',num2str(firstRow-2)];
Header4Range=['F',num2str(firstRow-2)];




writematrix(a1,xlsfilename,'Sheet',1,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',1,'Range',Header1Range)
writematrix('A1_LD_gp0',xlsfilename,'Sheet',1,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',1,'Range',Header3Range)
writematrix(char(Marker_address(1)),xlsfilename,'Sheet',1,'Range',Header4Range)


writematrix(a2,xlsfilename,'Sheet',2,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',2,'Range',Header1Range)
writematrix('A2_Rsquared_gp0',xlsfilename,'Sheet',2,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',2,'Range',Header3Range)
writematrix(char(Marker_address(2)),xlsfilename,'Sheet',2,'Range',Header4Range)


writematrix(a3,xlsfilename,'Sheet',3,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',3,'Range',Header1Range)
writematrix('A3_Reflectance_gp0',xlsfilename,'Sheet',3,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',3,'Range',Header3Range)
writematrix(char(Marker_address(3)),xlsfilename,'Sheet',3,'Range',Header4Range)




writematrix(a4,xlsfilename,'Sheet',4,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',4,'Range',Header1Range)
writematrix('A4_RMS_gp0',xlsfilename,'Sheet',4,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',4,'Range',Header3Range)
writematrix(char(Marker_address(4)),xlsfilename,'Sheet',4,'Range',Header4Range)



writematrix(a5,xlsfilename,'Sheet',5,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',5,'Range',Header1Range)
writematrix('A5_RmsNorm_gp0',xlsfilename,'Sheet',5,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',5,'Range',Header3Range)
writematrix(char(Marker_address(5)),xlsfilename,'Sheet',5,'Range',Header4Range)


writematrix(a6,xlsfilename,'Sheet',6,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',6,'Range',Header1Range)
writematrix('A6_RmsPoly_gp0',xlsfilename,'Sheet',6,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',6,'Range',Header3Range)
writematrix(char(Marker_address(6)),xlsfilename,'Sheet',6,'Range',Header4Range)


writematrix(a7,xlsfilename,'Sheet',7,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',7,'Range',Header1Range)
writematrix('A7_RmsSquared_gp0',xlsfilename,'Sheet',7,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',7,'Range',Header3Range)
writematrix(char(Marker_address(7)),xlsfilename,'Sheet',7,'Range',Header4Range)


writematrix(a8,xlsfilename,'Sheet',8,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',8,'Range',Header1Range)
writematrix('A8_Slope_gp0',xlsfilename,'Sheet',8,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',8,'Range',Header3Range)
writematrix(char(Marker_address(8)),xlsfilename,'Sheet',8,'Range',Header4Range)


writematrix(a9,xlsfilename,'Sheet',9,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',9,'Range',Header1Range)
writematrix('A9_LD_gp2',xlsfilename,'Sheet',9,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',9,'Range',Header3Range)
writematrix(char(Marker_address(9)),xlsfilename,'Sheet',9,'Range',Header4Range)


writematrix(a10,xlsfilename,'Sheet',10,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',10,'Range',Header1Range)
writematrix('A10_Rsquared_gp2',xlsfilename,'Sheet',10,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',10,'Range',Header3Range)
writematrix(char(Marker_address(10)),xlsfilename,'Sheet',10,'Range',Header4Range)


writematrix(a11,xlsfilename,'Sheet',11,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',11,'Range',Header1Range)
writematrix('A11_Reflectance_gp2',xlsfilename,'Sheet',11,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',11,'Range',Header3Range)
writematrix(char(Marker_address(11)),xlsfilename,'Sheet',11,'Range',Header4Range)


writematrix(a12,xlsfilename,'Sheet',12,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',12,'Range',Header1Range)
writematrix('A12_Rms_gp2',xlsfilename,'Sheet',12,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',12,'Range',Header3Range)
writematrix(char(Marker_address(12)),xlsfilename,'Sheet',12,'Range',Header4Range)


writematrix(a13,xlsfilename,'Sheet',13,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',13,'Range',Header1Range)
writematrix('A13_RmsNorm_gp2',xlsfilename,'Sheet',13,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',13,'Range',Header3Range)
writematrix(char(Marker_address(13)),xlsfilename,'Sheet',13,'Range',Header4Range)


writematrix(a14,xlsfilename,'Sheet',14,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',14,'Range',Header1Range)
writematrix('A14_RmsPoly_gp2',xlsfilename,'Sheet',14,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',14,'Range',Header3Range)
writematrix(char(Marker_address(14)),xlsfilename,'Sheet',14,'Range',Header4Range)



writematrix(a15,xlsfilename,'Sheet',15,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',15,'Range',Header1Range)
writematrix('A15_RmsSquared_gp2',xlsfilename,'Sheet',15,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',15,'Range',Header3Range)
writematrix(char(Marker_address(15)),xlsfilename,'Sheet',15,'Range',Header4Range)



writematrix(a16,xlsfilename,'Sheet',16,'Range',cellRange)
writematrix(MarkersName,xlsfilename,'Sheet',16,'Range',Header1Range)
writematrix('A16_Slope_gp2',xlsfilename,'Sheet',16,'Range',Header2Range)
writematrix(f1,xlsfilename,'Sheet',16,'Range',Header3Range)
writematrix(char(Marker_address(16)),xlsfilename,'Sheet',16,'Range',Header4Range)


end

