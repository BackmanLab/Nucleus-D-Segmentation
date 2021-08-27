function [MakerArray]=PWSLocalD100Markers(CellImage,NucMask)
clearvars -except NCOutput NCOutputCell NCOutputNuc ImageCell Mask NucMask CellImage
%Mask is the mask that select Nucleus boundry
%ImageCell is the image of Cell
Mask=NucMask;
ImageCell=CellImage;
%Read ROI MASK
NucImage=ImageCell;
NucImage(Mask==0)=0; %Creat an image of cell nucleus
%m=mean(mean(ImageCell(ImageCell>0)));
%s=std(mean(ImageCell));
mean_nuc=mean(NucImage(NucImage>0)); %Calculate the mean of the nucleus (the first threshold used in Fig1)
std_nuc=std(NucImage(NucImage>0));   %Calculate the standard deviation of the nucleus 
Marker43=mean_nuc;% Marker43 is Ave RMS of nuclues

Nonzero_Nuc=NucImage(NucImage>0);% all nucleus D values in 1D variable;
Sort_Nonzero_Nuc=sort(Nonzero_Nuc); %Sort the Nonzero value to find the 95% CI and the range for thresholding
Nuc_Length=max(size(Sort_Nonzero_Nuc));% Find the number of pixels in the Nuc. This Length help us to find the 95% CI of D in the Nucleus
Interval_25=round(Nuc_Length*2.5/100,0);% 2.5% of top and bottom D values which we use to find our 95% CI
LowD=Sort_Nonzero_Nuc(Interval_25);
HighD=Sort_Nonzero_Nuc(Nuc_Length-Interval_25);

%We determine four different thesohods  based on the mean of RMS in entire cell and standard deviation of RMS image in entire cell.  
ImageCell(ImageCell==min(min(ImageCell)))=0;
mean_cell1=mean(ImageCell(ImageCell>0));
%mean_cell1=mean(ImageCell(ImageCell>0)); %We calculate the first threshold based on the mean of entire cell.
%Threshold_High=std(ImageCell(ImageCell>0));    %We calculate the gap between different threshold levels based on the standard deviation of entire cell image. 
Threshold_High= (HighD-mean_cell1)/5;
Threshold_Low= (mean_cell1 - LowD)/5;


mean_cell2=mean_cell1+Threshold_High;         %
mean_cell3=mean_cell1+2*Threshold_High;
mean_cell4=mean_cell1+3*Threshold_High;
mean_cell5=mean_cell1+4*Threshold_High;

mean_cell6=mean_cell1 - Threshold_High;
mean_cell7=mean_cell1 - 2*Threshold_High;



NucArea=NucImage;
NucArea(NucImage>0)=1;%highlight the Area of the nucleus with 1 and other area with zero 
%%%subplot(3,3,2);
%%%imshow(NucImage)% Plot only Nucleus. The firs image we show is the nucleus
%%%fontsize=16;
%%%Margin1=76;
%%%Margin2=100;
%%%Margin3=200;

%%%text_str = ['Nucleus Ave RMS: ' num2str(mean_nuc,'%0.2f') '%'];
%%%hText = text(Margin1,Margin2,text_str,'Color',[1 1 0],'FontSize',fontsize);

%--------------------------
%Biomarkers1,2,3,4,5,6,
% We define six biomarkers.

%Biomarker1 is "HighRMSMean1" which is the
% average RMS value of areas of nuclues that have an RMS greater than
% the "average RMS of entire cell".

%Biomarker2 is "HighRMSAreaRatio1" which
% is the ratio of area of the nuclues with average RMS > average RMS of the
% entire cell to the total area of the nuclues.

% Biomarker3 is "HighRMSSTD1" which is Standard deviation between all the pixels greater than "mean_cell1" Threshold.

%Biomarker4 is "HighRMSMean2" which is the Average RMS in the area of the cells with RMS value winthin the following range of "mean_cell1" <RMS< "mean_cell2" 

%Biomarker5 is "HighRMSAreaRatio2" which is the is the ratio of area of the nuclues with average RMS greater average RMS of the
% entire cell but smaller than average RMS of entire cell +one standard deviation to the total area of the nuclues.

% Biomarker6 is "HighRMSSTD2" which is std of RMS in the locations withaverage RMS greater average RMS of the
% entire cell but smaller than average RMS of entire cell +one standard deviation 



%%%subplot(3,3,4);
NucHigh1=NucImage;
NucHigh1(NucHigh1<mean_cell1)=0;%Biomarker formation
%%%imshow(NucHigh2)
NucHigh1Area=NucHigh1;
NucHigh1Area(NucHigh1Area>0)=1;%Biomarker formation

NucHigh2=NucHigh1;
NucHigh2(NucHigh2>mean_cell2)=0;% Biomarker formation...Nuclues area greater than Mean_Cell1 and smaller than Mean_Cell2

NucHigh2Area=NucHigh2;
NucHigh2Area(NucHigh2Area>0)=1;%Biomarker formation


HighRMSMean1= sum(sum(NucHigh1))/sum(sum(NucHigh1Area));%BioMarker1:Average RMS in the High RMS Area (determined by "mean_cell1" Threshold)
HighRMSAreaRatio1=100*(sum(sum(NucHigh1Area))/sum(sum(NucArea)));%BioMarker2: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" Threshold)
HighRMSSTD1=std(NucHigh1(NucHigh1>0));% Biomarker3 is Standard deviation between all the pixels greater than "mean_cell1" Threshold.

HighRMSMean2= sum(sum(NucHigh2))/sum(sum(NucHigh2Area));%BioMarker4:Average RMS in the High RMS Area (determined by "mean_cell1" <RMS< "mean_cell2" Threshold)
HighRMSAreaRatio2=100*(sum(sum(NucHigh2Area))/sum(sum(NucArea)));%BioMarker5: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" <RMS< "mean_cell2" Threshold)
HighRMSSTD2=std(NucHigh2(NucHigh2>0)); %BioMarker6: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" <RMS< "mean_cell2" Threshold)


%%%text_str = ['Nucleus High RMS2 Mean: ' num2str(HighRMSMean2,'%0.2f') '%'];
%%%hText = text(Margin1,Margin2,text_str,'Color',[1 1 0],'FontSize',fontsize);
%%%text_str = ['Nucleus HighRMS2 Area: ' num2str(HighRMSAreaRatio2,'%0.2f') '%'];
%%%hText = text(Margin1,Margin3,text_str,'Color',[1 1 0],'FontSize',fontsize);
%----------------------------------------------
%Biomarkers7,8,9,10,11,12,
% We define six biomarkers.

%Biomarker7 is "HighRMSMean3" which is the
% average RMS value of areas of nuclues that have an RMS greater than
% the "average RMS of entire cell +one STD".

%Biomarker8 is "HighRMSAreaRatio3" which
% is the ratio of area of the nuclues with average RMS > average RMS of the
% entire cell +one STD to the total area of the nuclues.

% Biomarker9 is "HighRMSSTD3" which is Standard deviation between all the pixels greater than "mean_cell2" Threshold.

%Biomarker10 is "HighRMSMean4" which is the Average RMS in the area of the cells with RMS value winthin the following range of "mean_cell2" <RMS< "mean_cell3" 

%Biomarker11 is "HighRMSAreaRatio4" which is the is the ratio of area of the nuclues with average RMS greater than the "average RMS of the
% entire cell + one STD" but smaller than the "average RMS of entire cell + two standard deviation to the total area of the nuclues".

% Biomarker12 is "HighRMSSTD4" which is std of RMS in the locations withaverage RMS greater average RMS of the
% entire cell + one STD but smaller than average RMS of entire cell + two standard deviation 




%%%subplot(3,3,4);
NucHigh3=NucImage;
NucHigh3(NucHigh3<mean_cell2)=0;% Biomarker formation...For Biomarker7
%%%imshow(NucHigh2)
NucHigh3Area=NucHigh3;
NucHigh3Area(NucHigh3Area>0)=1;% Biomarker formation...For Biomarker8

NucHigh4=NucHigh3;
NucHigh4(NucHigh4>mean_cell3)=0;% Biomarker formation...For Biomarker9. Nuclues area greater than Mean_Cell2 and smaller than Mean_Cell3

NucHigh4Area=NucHigh4;
NucHigh4Area(NucHigh4Area>0)=1;%Biomarker formation...For Biomarker10.


HighRMSMean3= sum(sum(NucHigh3))/sum(sum(NucHigh3Area));%BioMarker7:Average RMS in the High RMS Area (determined by "mean_cell1" Threshold)
HighRMSAreaRatio3=100*(sum(sum(NucHigh3Area))/sum(sum(NucArea)));%BioMarker8: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" Threshold)
HighRMSSTD3=std(NucHigh3(NucHigh3>0));%Biomarker 9

HighRMSMean4= sum(sum(NucHigh4))/sum(sum(NucHigh4Area));%BioMarker10:Average RMS in the High RMS Area (determined by "mean_cell1" <RMS< "mean_cell2" Threshold)
HighRMSAreaRatio4=100*(sum(sum(NucHigh4Area))/sum(sum(NucArea)));%BioMarker11: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" <RMS< "mean_cell2" Threshold)
HighRMSSTD4=std(NucHigh4(NucHigh4>0));%Biomarker12






% End for Biomarkers 7,8,9,10,11,12
%----------------------------------------------

%----------------------------------------------
%Biomarkers13,14,15,16,17,18,
% We define six biomarkers.

%Biomarker13 is "HighRMSMean5" which is the
% average RMS value of areas of nuclues that have an RMS greater than
% the "average RMS of entire cell + two STD".

%Biomarker14 is "HighRMSAreaRatio5" which
% is the ratio of area of the nuclues with average RMS > average RMS of the
% entire cell +two STD to the total area of the nuclues.

% Biomarker15 is "HighRMSSTD5" which is Standard deviation between all the pixels greater than "mean_cell3" Threshold.

%Biomarker16 is "HighRMSMean6" which is the Average RMS in the area of the cells with RMS value winthin the following range of "mean_cell3" <RMS< "mean_cell4" 

%Biomarker17 is "HighRMSAreaRatio6" which is the is the ratio of area of the nuclues with average RMS greater than the "average RMS of the
% entire cell + two STD" but smaller than the "average RMS of entire cell + three standard deviation to the total area of the nuclues".

% Biomarker18 is "HighRMSSTD6" which is std of RMS in the locations withaverage RMS greater average RMS of the
% entire cell + two STD but smaller than average RMS of entire cell + three standard deviation 



%%%subplot(3,3,4);
NucHigh5=NucImage;
NucHigh5(NucHigh5<mean_cell3)=0;% For Biomarker9
%%%imshow(NucHigh2)
NucHigh5Area=NucHigh5;
NucHigh5Area(NucHigh5Area>0)=1;% For Biomarker10


NucHigh6=NucHigh5;
NucHigh6(NucHigh6>mean_cell4)=0;% For Biomarker11. Nuclues area greater than Mean_Cell2 and smaller than Mean_Cell3

NucHigh6Area=NucHigh6;
NucHigh6Area(NucHigh6Area>0)=1;%For Biomarker12.


HighRMSMean5= sum(sum(NucHigh5))/sum(sum(NucHigh5Area));%BioMarker1:Average RMS in the High RMS Area (determined by "mean_cell1" Threshold)
HighRMSAreaRatio5=100*(sum(sum(NucHigh5Area))/sum(sum(NucArea)));%BioMarker2: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" Threshold)
HighRMSSTD5=std(NucHigh5(NucHigh5>0));

HighRMSMean6= sum(sum(NucHigh6))/sum(sum(NucHigh6Area));%BioMarker3:Average RMS in the High RMS Area (determined by "mean_cell1" <RMS< "mean_cell2" Threshold)
HighRMSAreaRatio6=100*(sum(sum(NucHigh6Area))/sum(sum(NucArea)));%BioMarker4: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" <RMS< "mean_cell2" Threshold)
HighRMSSTD6=std(NucHigh6(NucHigh6>0));

% End for Biomarkers 13,14,15,16,17,18,
%----------------------------------------------


%----------------------------------------------
%Biomarkers19,20,21,22,23,24,
% We define six biomarkers.

%Biomarker19 is "HighRMSMean7" which is the
% average RMS value of areas of nuclues that have an RMS greater than
% the "average RMS of entire cell + three STD".

%Biomarker20 is "HighRMSAreaRatio7" which
% is the ratio of area of the nuclues with average RMS > average RMS of the
% entire cell +three STD to the total area of the nuclues.

% Biomarker21 is "HighRMSSTD7" which is Standard deviation between all the pixels greater than "mean_cell4" Threshold.

%Biomarker22 is "HighRMSMean8" which is the Average RMS in the area of the cells with RMS value winthin the following range of "mean_cell4" <RMS< "mean_cell5" 

%Biomarker23 is "HighRMSAreaRatio8" which is the is the ratio of area of the nuclues with average RMS greater than the "average RMS of the
% entire cell + three STD" but smaller than the "average RMS of entire cell + four standard deviation to the total area of the nuclues".

% Biomarker24 is "HighRMSSTD8" which is std of RMS in the locations withaverage RMS greater average RMS of the
% entire cell + three STD but smaller than average RMS of entire cell + four standard deviation 



%%%subplot(3,3,4);
NucHigh7=NucImage;
NucHigh7(NucHigh7<mean_cell4)=0;% For Biomarker13
%%%imshow(NucHigh2)
NucHigh7Area=NucHigh7;
NucHigh7Area(NucHigh7Area>0)=1;% For Biomarker14

NucHigh8=NucHigh7;
NucHigh8(NucHigh8>mean_cell5)=0;% For Biomarker15. Nuclues area greater than Mean_Cell2 and smaller than Mean_Cell3

NucHigh8Area=NucHigh8;
NucHigh8Area(NucHigh8Area>0)=1;%For Biomarker16.


HighRMSMean7= sum(sum(NucHigh7))/sum(sum(NucHigh7Area));%BioMarker1:Average RMS in the High RMS Area (determined by "mean_cell1" Threshold)
HighRMSAreaRatio7=100*(sum(sum(NucHigh7Area))/sum(sum(NucArea)));%BioMarker2: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" Threshold)
HighRMSSTD7=std(NucHigh7(NucHigh7>0));

HighRMSMean8= sum(sum(NucHigh8))/sum(sum(NucHigh8Area));%BioMarker3:Average RMS in the High RMS Area (determined by "mean_cell1" <RMS< "mean_cell2" Threshold)
HighRMSAreaRatio8=100*(sum(sum(NucHigh8Area))/sum(sum(NucArea)));%BioMarker4: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" <RMS< "mean_cell2" Threshold)
HighRMSSTD8=std(NucHigh8(NucHigh8>0));
% End for Biomarkers 19,20,21,22,23,24,

%-------------------------

%----------------------------------------------
%Biomarkers25,26,27,28,29,30
% We define six biomarkers.

%Biomarker25 is "HighRMSMean9" which is the
% average RMS value of areas of nuclues that have an RMS greater than
% the "average RMS of entire cell + three STD".

%Biomarker26 is "HighRMSAreaRatio9" which
% is the ratio of area of the nuclues with average RMS > average RMS of the
% entire cell +four STD to the total area of the nuclues.

% Biomarker27 is "HighRMSSTD9" which is Standard deviation between all the pixels greater than "mean_cell5" Threshold.

%Biomarker28 is "HighRMSMean10" which is the Average RMS of area of nuclues
%with where each pixel is greater than nuclues average RMS

%Biomarker29 is "HighRMSAreaRatio10" whic his the ratio of area of the nuclues with average RMS > average RMS of the
% entire nuclues

%Biomarker30 is "HighRMSSTD10" is Standard deviation between all the pixels greater with average RMS > average RMS of the
% entire nuclues

%%%subplot(3,3,4);
NucHigh9=NucImage;
NucHigh9(NucHigh9<mean_cell5)=0;% For Biomarker17
%%%imshow(NucHigh2)
NucHigh9Area=NucHigh9;
NucHigh9Area(NucHigh9Area>0)=1;% For Biomarker18



HighRMSMean9= sum(sum(NucHigh9))/sum(sum(NucHigh9Area));%BioMarker1:Average RMS in the High RMS Area (determined by "mean_cell1" Threshold)
HighRMSAreaRatio9=100*(sum(sum(NucHigh9Area))/sum(sum(NucArea)));%BioMarker2: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" Threshold)
HighRMSSTD9=std(NucHigh9(NucHigh9>0));

% End for Biomarkers 17,18

%-----------------------------------------------

%-----------------------------
%We only select the are of nucleus that have an RMS>average nucleus RMS
%%%subplot(3,3,3);
NucHigh10=NucImage;% NucHigh1 become the RMS image of nuclues
NucHigh10(NucHigh10<mean_nuc)=0;% NucHigh1 becomes the RMS image of nuclues that is greater average of nuclues RMS.
%%%imshow(NucHigh1)



NucHigh10Area=NucHigh10;
NucHigh10Area(NucHigh10Area>0)=1;% Calculate the area of nuclues that has RMS > average Nuclues RMS

%Calculate average RMS of HighRMS AREA using threshold 1 = mean of the nucleus RMS 
HighRMSMean10= sum(sum(NucHigh10))/sum(sum(NucHigh10Area));%Marker 1: Average RMS of High RMS Area of nucleus(Threshold 1 = mean of the nucleus)
%Calculate ratio of HighRMS AREA to LowRMS ares using threshold 1=mean of the nucleus RMS 
HighRMSAreaRatio10=100*(sum(sum(NucHigh10Area))/sum(sum(NucArea))); % Marker 2: HIGH RMS Area devided b Low RMS Area (Threshold 1 = mean of the nucleus)
HighRMSSTD10=std(NucHigh10(NucHigh10>0));

%%%text_str = ['Nucleus High RMS1 Mean: ' num2str(HighRMSMean1,'%0.2f') '%'];
%%%hText = text(Margin1,Margin2,text_str,'Color',[1 1 0],'FontSize',fontsize);

%%%text_str = ['Nucleus HighRMS1 Area: ' num2str(HighRMSAreaRatio1,'%0.2f') '%'];
%%%hText = text(Margin1,200,text_str,'Color',[1 1 0],'FontSize',fontsize);



%--------------- Mean (Cell_Rms) - STD < Rms < Mean (Cell_Rms)
%--------------Start
%%%subplot(3,3,4);
NucHigh11=NucImage;
NucHigh11(NucHigh11>mean_cell1)=0;%Biomarker formation
%%%imshow(NucHigh2)
NucHigh11Area=NucHigh11;
NucHigh11Area(NucHigh11Area>0)=1;%Biomarker formation

NucHigh12=NucHigh11;
NucHigh12(NucHigh12 < mean_cell6)=0;% Biomarker formation...Nuclues area greater than Mean_Cell1 and smaller than Mean_Cell2

NucHigh12Area=NucHigh12;
NucHigh12Area(NucHigh12Area>0)=1;%Biomarker formation


HighRMSMean11= sum(sum(NucHigh11))/sum(sum(NucHigh11Area));%BioMarker1:Average RMS in the High RMS Area (determined by "mean_cell1" Threshold)
HighRMSAreaRatio11=100*(sum(sum(NucHigh11Area))/sum(sum(NucArea)));%BioMarker2: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" Threshold)
HighRMSSTD11=std(NucHigh11(NucHigh11>0));% Biomarker3 is Standard deviation between all the pixels greater than "mean_cell1" Threshold.

HighRMSMean12= sum(sum(NucHigh12))/sum(sum(NucHigh12Area));%BioMarker4:Average RMS in the High RMS Area (determined by "mean_cell1" <RMS< "mean_cell2" Threshold)
HighRMSAreaRatio12=100*(sum(sum(NucHigh12Area))/sum(sum(NucArea)));%BioMarker5: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" <RMS< "mean_cell2" Threshold)
HighRMSSTD12=std(NucHigh12(NucHigh12>0)); %BioMarker6: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" <RMS< "mean_cell2" Threshold)


%%%text_str = ['Nucleus High RMS2 Mean: ' num2str(HighRMSMean2,'%0.2f') '%'];
%%%hText = text(Margin1,Margin2,text_str,'Color',[1 1 0],'FontSize',fontsize);
%%%text_str = ['Nucleus HighRMS2 Area: ' num2str(HighRMSAreaRatio2,'%0.2f') '%'];
%%%hText = text(Margin1,Margin3,text_str,'Color',[1 1 0],'FontSize',fontsize);
%----------------------------------------------
%Biomarkers7,8,9,10,11,12,
% We define six biomarkers.

%Biomarker7 is "HighRMSMean3" which is the
% average RMS value of areas of nuclues that have an RMS greater than
% the "average RMS of entire cell +one STD".

%Biomarker8 is "HighRMSAreaRatio3" which
% is the ratio of area of the nuclues with average RMS > average RMS of the
% entire cell +one STD to the total area of the nuclues.

% Biomarker9 is "HighRMSSTD3" which is Standard deviation between all the pixels greater than "mean_cell2" Threshold.

%Biomarker10 is "HighRMSMean4" which is the Average RMS in the area of the cells with RMS value winthin the following range of "mean_cell2" <RMS< "mean_cell3" 

%Biomarker11 is "HighRMSAreaRatio4" which is the is the ratio of area of the nuclues with average RMS greater than the "average RMS of the
% entire cell + one STD" but smaller than the "average RMS of entire cell + two standard deviation to the total area of the nuclues".

% Biomarker12 is "HighRMSSTD4" which is std of RMS in the locations withaverage RMS greater average RMS of the
% entire cell + one STD but smaller than average RMS of entire cell + two standard deviation 




%%%subplot(3,3,4);
NucHigh13=NucImage;
NucHigh13(NucHigh13>mean_cell6)=0;% Biomarker formation...For Biomarker7
%%%imshow(NucHigh2)
NucHigh13Area=NucHigh13;
NucHigh13Area(NucHigh13Area>0)=1;% Biomarker formation...For Biomarker8

NucHigh14=NucHigh13;
NucHigh14(NucHigh14<mean_cell7)=0;% Biomarker formation...For Biomarker9. Nuclues area greater than Mean_Cell2 and smaller than Mean_Cell3

NucHigh14Area=NucHigh14;
NucHigh14Area(NucHigh14Area>0)=1;%Biomarker formation...For Biomarker10.


HighRMSMean13= sum(sum(NucHigh13))/sum(sum(NucHigh13Area));%BioMarker1:Average RMS in the High RMS Area (determined by "mean_cell1" Threshold)
HighRMSAreaRatio13=100*(sum(sum(NucHigh13Area))/sum(sum(NucArea)));%BioMarker2: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" Threshold)
HighRMSSTD13=std(NucHigh13(NucHigh13>0));% Biomarker3 is Standard deviation between all the pixels greater than "mean_cell1" Threshold.

HighRMSMean14= sum(sum(NucHigh14))/sum(sum(NucHigh14Area));%BioMarker4:Average RMS in the High RMS Area (determined by "mean_cell1" <RMS< "mean_cell2" Threshold)
HighRMSAreaRatio14=100*(sum(sum(NucHigh14Area))/sum(sum(NucArea)));%BioMarker5: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" <RMS< "mean_cell2" Threshold)
HighRMSSTD14=std(NucHigh14(NucHigh14>0)); %BioMarker6: Area of High RMS devided by Nucleus Area(determined by "mean_cell1" <RMS< "mean_cell2" Threshold)


NucArea1Connected1= bwpropfilt(logical(NucHigh1Area),'perimeter',1);
NucArea1Connected2= bwpropfilt(logical(NucHigh1Area),'perimeter',2);
NucArea1Connected3= bwpropfilt(logical(NucHigh1Area),'perimeter',3);

NucArea2Connected1= bwpropfilt(logical(NucHigh2Area),'perimeter',1);
NucArea2Connected2= bwpropfilt(logical(NucHigh2Area),'perimeter',2);
NucArea2Connected3= bwpropfilt(logical(NucHigh2Area),'perimeter',3);

NucArea3Connected1= bwpropfilt(logical(NucHigh3Area),'perimeter',1);
NucArea3Connected2= bwpropfilt(logical(NucHigh3Area),'perimeter',2);
NucArea3Connected3= bwpropfilt(logical(NucHigh3Area),'perimeter',3);

NucArea4Connected1= bwpropfilt(logical(NucHigh4Area),'perimeter',1);
NucArea4Connected2= bwpropfilt(logical(NucHigh4Area),'perimeter',2);
NucArea4Connected3= bwpropfilt(logical(NucHigh4Area),'perimeter',3);

NucArea5Connected1= bwpropfilt(logical(NucHigh5Area),'perimeter',1);
NucArea5Connected2= bwpropfilt(logical(NucHigh5Area),'perimeter',2);
NucArea5Connected3= bwpropfilt(logical(NucHigh5Area),'perimeter',3);

NucArea6Connected1= bwpropfilt(logical(NucHigh6Area),'perimeter',1);
NucArea6Connected2= bwpropfilt(logical(NucHigh6Area),'perimeter',2);
NucArea6Connected3= bwpropfilt(logical(NucHigh6Area),'perimeter',3);

NucArea7Connected1= bwpropfilt(logical(NucHigh7Area),'perimeter',1);
NucArea7Connected2= bwpropfilt(logical(NucHigh7Area),'perimeter',2);
NucArea7Connected3= bwpropfilt(logical(NucHigh7Area),'perimeter',3);


NucArea8Connected1= bwpropfilt(logical(NucHigh8Area),'perimeter',1);
NucArea8Connected2= bwpropfilt(logical(NucHigh8Area),'perimeter',2);
NucArea8Connected3= bwpropfilt(logical(NucHigh8Area),'perimeter',3);


NucArea9Connected1= bwpropfilt(logical(NucHigh9Area),'perimeter',1);
NucArea9Connected2= bwpropfilt(logical(NucHigh9Area),'perimeter',2);
NucArea9Connected3= bwpropfilt(logical(NucHigh9Area),'perimeter',3);


NucArea11Connected1= bwpropfilt(logical(NucHigh11Area),'perimeter',1);
NucArea11Connected2= bwpropfilt(logical(NucHigh11Area),'perimeter',2);
NucArea11Connected3= bwpropfilt(logical(NucHigh11Area),'perimeter',3);


NucArea12Connected1= bwpropfilt(logical(NucHigh12Area),'perimeter',1);
NucArea12Connected2= bwpropfilt(logical(NucHigh12Area),'perimeter',2);
NucArea12Connected3= bwpropfilt(logical(NucHigh12Area),'perimeter',3);


NucArea13Connected1= bwpropfilt(logical(NucHigh13Area),'perimeter',1);
NucArea13Connected2= bwpropfilt(logical(NucHigh13Area),'perimeter',2);
NucArea13Connected3= bwpropfilt(logical(NucHigh13Area),'perimeter',3);


NucArea14Connected1= bwpropfilt(logical(NucHigh14Area),'perimeter',1);
NucArea14Connected2= bwpropfilt(logical(NucHigh14Area),'perimeter',2);
NucArea14Connected3= bwpropfilt(logical(NucHigh14Area),'perimeter',3);

%----------------End

%%%%Boundy and Center Related Markers

[Marker_Boundry1 Marker_Boundry2 Marker_Center]= BoundryMarkers(NucImage,Mask,5);



%End
%---------------------------------
Marker1= HighRMSMean1;
Marker2= HighRMSAreaRatio1;
Marker3= HighRMSSTD1;

Marker4= HighRMSMean2;
Marker5= HighRMSAreaRatio2;
Marker6= HighRMSSTD2;

Marker7= HighRMSMean3;
Marker8= HighRMSAreaRatio3;
Marker9= HighRMSSTD3;

Marker10= HighRMSMean4;
Marker11= HighRMSAreaRatio4;
Marker12= HighRMSSTD4;

Marker13= HighRMSMean5;
Marker14= HighRMSAreaRatio5;
Marker15= HighRMSSTD5;

Marker16= HighRMSMean6;
Marker17= HighRMSAreaRatio6;
Marker18= HighRMSSTD6;

Marker19= HighRMSMean7;
Marker20= HighRMSAreaRatio7;
Marker21= HighRMSSTD7;

Marker22= HighRMSMean8;
Marker23= HighRMSAreaRatio8;
Marker24= HighRMSSTD8;

Marker25= HighRMSMean9;
Marker26= HighRMSAreaRatio9;
Marker27= HighRMSSTD9;

Marker28= HighRMSMean10;
Marker29= HighRMSAreaRatio10;
Marker30= HighRMSSTD10;


Marker31=HighRMSMean3-HighRMSMean1;
Marker32=HighRMSAreaRatio1-HighRMSAreaRatio3;

Marker33=HighRMSMean5-HighRMSMean3;
Marker34=HighRMSAreaRatio3-HighRMSAreaRatio5;

Marker35=HighRMSMean7-HighRMSMean5;
Marker36=HighRMSAreaRatio5-HighRMSAreaRatio7;

Marker37=HighRMSMean9-HighRMSMean7;
Marker38=HighRMSAreaRatio7-HighRMSAreaRatio9;

Marker39=std([Marker31 Marker33 Marker35 Marker37]);
Marker40=std([Marker32 Marker34 Marker36 Marker38]);

Marker41=(abs(Marker31)+abs(Marker33)+abs(Marker35)+abs(Marker37))/4;
Marker42=(abs(Marker32)+abs(Marker34)+abs(Marker36)+abs(Marker38))/4;
Marker43=Marker43;

Marker44= HighRMSMean11;
Marker45= HighRMSAreaRatio11;
Marker46= HighRMSSTD11;

Marker47= HighRMSMean12;
Marker48= HighRMSAreaRatio12;
Marker49= HighRMSSTD12;

Marker50= HighRMSMean13;
Marker51= HighRMSAreaRatio13;
Marker52= HighRMSSTD13;

Marker53= HighRMSMean14;
Marker54= HighRMSAreaRatio14;
Marker55= HighRMSSTD14;

%Maybe and Mayb not

Marker56=100*(sum(sum(double(NucArea1Connected1)))/sum(sum(NucArea)));
Marker57=100*(sum(sum(double(NucArea1Connected2)))/sum(sum(NucArea)));
Marker58=100*(sum(sum(double(NucArea1Connected3)))/sum(sum(NucArea)));

Marker59=100*(sum(sum(double(NucArea2Connected1)))/sum(sum(NucArea)));
Marker60=100*(sum(sum(double(NucArea2Connected2)))/sum(sum(NucArea)));
Marker61=100*(sum(sum(double(NucArea2Connected3)))/sum(sum(NucArea)));

Marker62=100*(sum(sum(double(NucArea3Connected1)))/sum(sum(NucArea)));
Marker63=100*(sum(sum(double(NucArea3Connected2)))/sum(sum(NucArea)));
Marker64=100*(sum(sum(double(NucArea3Connected3)))/sum(sum(NucArea)));

Marker65=100*(sum(sum(double(NucArea4Connected1)))/sum(sum(NucArea)));
Marker66=100*(sum(sum(double(NucArea4Connected2)))/sum(sum(NucArea)));
Marker67=100*(sum(sum(double(NucArea4Connected3)))/sum(sum(NucArea)));

Marker68=100*(sum(sum(double(NucArea5Connected1)))/sum(sum(NucArea)));
Marker69=100*(sum(sum(double(NucArea5Connected2)))/sum(sum(NucArea)));
Marker70=100*(sum(sum(double(NucArea5Connected3)))/sum(sum(NucArea)));

Marker71=100*(sum(sum(double(NucArea6Connected1)))/sum(sum(NucArea)));
Marker72=100*(sum(sum(double(NucArea6Connected2)))/sum(sum(NucArea)));
Marker73=100*(sum(sum(double(NucArea6Connected3)))/sum(sum(NucArea)));

Marker74=100*(sum(sum(double(NucArea7Connected1)))/sum(sum(NucArea)));
Marker75=100*(sum(sum(double(NucArea7Connected2)))/sum(sum(NucArea)));
Marker76=100*(sum(sum(double(NucArea7Connected3)))/sum(sum(NucArea)));

Marker77=100*(sum(sum(double(NucArea8Connected1)))/sum(sum(NucArea)));
Marker78=100*(sum(sum(double(NucArea8Connected2)))/sum(sum(NucArea)));
Marker79=100*(sum(sum(double(NucArea8Connected3)))/sum(sum(NucArea)));

Marker80=100*(sum(sum(double(NucArea9Connected1)))/sum(sum(NucArea)));
Marker81=100*(sum(sum(double(NucArea9Connected2)))/sum(sum(NucArea)));
Marker82=100*(sum(sum(double(NucArea9Connected3)))/sum(sum(NucArea)));

Marker83=100*(sum(sum(double(NucArea11Connected1)))/sum(sum(NucArea)));
Marker84=100*(sum(sum(double(NucArea11Connected2)))/sum(sum(NucArea)));
Marker85=100*(sum(sum(double(NucArea11Connected3)))/sum(sum(NucArea)));

Marker86=100*(sum(sum(double(NucArea12Connected1)))/sum(sum(NucArea)));
Marker87=100*(sum(sum(double(NucArea12Connected2)))/sum(sum(NucArea)));
Marker88=100*(sum(sum(double(NucArea12Connected3)))/sum(sum(NucArea)));

Marker89=100*(sum(sum(double(NucArea13Connected1)))/sum(sum(NucArea)));
Marker90=100*(sum(sum(double(NucArea13Connected2)))/sum(sum(NucArea)));
Marker91=100*(sum(sum(double(NucArea13Connected3)))/sum(sum(NucArea)));

Marker92=100*(sum(sum(double(NucArea14Connected1)))/sum(sum(NucArea)));
Marker93=100*(sum(sum(double(NucArea14Connected2)))/sum(sum(NucArea)));
Marker94=100*(sum(sum(double(NucArea14Connected3)))/sum(sum(NucArea)));

Marker95=Marker_Boundry1;
Marker96=Marker_Boundry2;
Marker97=Marker_Center;
Marker98=Marker_Boundry1/Marker_Boundry2;
Marker99=Marker_Boundry1/Marker_Center;
Marker100=Marker_Boundry2/Marker_Center;



MakerArray=[Marker1 Marker2 Marker3 Marker4 Marker5 Marker6 Marker7 Marker8 Marker9 Marker10 Marker11 Marker12 Marker13 Marker14 Marker15 Marker16 Marker17 Marker18 Marker19 Marker20 Marker21 Marker22 Marker23 Marker24 Marker25 Marker26 Marker27 Marker28 Marker29 Marker30 Marker31 Marker32 Marker33 Marker34 Marker35 Marker36 Marker37 Marker38 Marker39 Marker40 Marker41 Marker42 Marker43 Marker44 Marker45 Marker46 Marker47 Marker48 Marker49 Marker50 Marker51 Marker52 Marker53 Marker54 Marker55 Marker56 Marker57 Marker58 Marker59 Marker60 Marker61 Marker62 Marker63 Marker64 Marker65 Marker66 Marker67 Marker68 Marker69 Marker70 Marker71 Marker72 Marker73 Marker74 Marker75 Marker76 Marker77 Marker78 Marker79 Marker80 Marker81 Marker82 Marker83 Marker84 Marker85 Marker86 Marker87 Marker88 Marker89 Marker90 Marker91 Marker92 Marker93 Marker94 Marker95 Marker96 Marker97 Marker98 Marker99 Marker100];

%Questions can be adressed to alid@northwestern.edu
