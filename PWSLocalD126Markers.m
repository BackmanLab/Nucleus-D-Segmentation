function [MakerArray]=PWSLocalD126Markers(CellImage,NucMask)
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
TraditionalRms=mean_nuc;% Marker43 is Ave RMS of nuclues

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
Threshold_Low= (mean_cell1 - LowD)/2;


mean_cell2=mean_cell1+Threshold_High;         %
mean_cell3=mean_cell1+2*Threshold_High;
mean_cell4=mean_cell1+3*Threshold_High;
mean_cell5=mean_cell1+4*Threshold_High;

mean_cell6=mean_cell1 - Threshold_Low;
mean_cell7=mean_cell1 - 2*Threshold_Low;



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
%Low- Packing domain nsubsections (we will present them in blue)
%%%subplot(3,3,3);
NucHigh10=NucImage;% NucHigh10 become the D image of nuclues
NucHigh10(NucHigh10>mean_cell7)=0;% NucHigh10 becomes the D-image subsection of nuclues that is smaller thanmean_cell7  average of nuclues RMS.
%%%imshow(NucHigh10)



NucHigh10Area=NucHigh10;
NucHigh10Area(NucHigh10Area>0)=1;% Calculate the area of nuclues subsection that D < MeanCell Image-2 threshold (mean_cell7)

HighRMSMean10= sum(sum(NucHigh10))/sum(sum(NucHigh10Area));%Marker 1: Average D value of nucleus subsection: D < MeanCell Image-2 threshold (mean_cell7)
HighRMSAreaRatio10=100*(sum(sum(NucHigh10Area))/sum(sum(NucArea))); % Marker 2: Area ratio of lowD subsection to the entire nucleus: D < MeanCell Image-2 threshold (mean_cell7)
HighRMSSTD10=std(NucHigh10(NucHigh10>0));% Marker 3: Standard deviation of lowD subsection to the entire nucleus: D < MeanCell Image-2 threshold (mean_cell7)

%%%text_str = ['Nucleus High RMS1 Mean: ' num2str(HighRMSMean1,'%0.2f') '%'];
%%%hText = text(Margin1,Margin2,text_str,'Color',[1 1 0],'FontSize',fontsize);

%%%text_str = ['Nucleus HighRMS1 Area: ' num2str(HighRMSAreaRatio1,'%0.2f') '%'];
%%%hText = text(Margin1,200,text_str,'Color',[1 1 0],'FontSize',fontsize);


%%%subplot(3,3,4);
NucHigh11=NucImage;
NucHigh11(NucHigh11>mean_cell6)=0;% For Marker 4: Average D value of nucleus subsection: D < MeanCell Image- threshold (mean_cell6)
%%%imshow(NucHigh2)
NucHigh11Area=NucHigh11;
NucHigh11Area(NucHigh11Area>0)=1;% % For Marker 5: Area ratio of lowD subsection to the entire nucleus: D < MeanCell Image- threshold (mean_cell)

NucHigh12=NucHigh11;
NucHigh12(NucHigh12<=mean_cell7)=0;% For Marker 7: Average D value of lowD nucleus subsection: MeanCell Image- 2*threshold (mean_cell7)< D < MeanCell Image- threshold (mean_cell6)

NucHigh12Area=NucHigh12;
NucHigh12Area(NucHigh12Area>0)=1;%For Marker 8: Area ratio of nucleus LowD subsection: MeanCell Image- 2*threshold (mean_cell7)< D < MeanCell Image- threshold (mean_cell6)


HighRMSMean11= sum(sum(NucHigh11))/sum(sum(NucHigh11Area));%Marker 4: Average D value of nucleus subsection: D < MeanCell Image- threshold (mean_cell6)
HighRMSAreaRatio11=100*(sum(sum(NucHigh11Area))/sum(sum(NucArea))); %Marker 5: Area ratio of lowD subsection to the entire nucleus: D < MeanCell Image- threshold (mean_cell)

HighRMSSTD11=std(NucHigh11(NucHigh11>0));% Marker 6: Nonzero-Satandard deviation of lowD subsection to the entire nucleus: D < MeanCell Image- threshold (mean_cell)

HighRMSMean12= sum(sum(NucHigh12))/sum(sum(NucHigh12Area));%Marker 7: Average D value of lowD nucleus subsection: MeanCell Image- 2*threshold (mean_cell7)< D < MeanCell Image- threshold (mean_cell6)
HighRMSAreaRatio12=100*(sum(sum(NucHigh12Area))/sum(sum(NucArea)));%Marker 8: Area ratio of nucleus LowD subsection: MeanCell Image- 2*threshold (mean_cell7)< D < MeanCell Image- threshold (mean_cell6)
HighRMSSTD12=std(NucHigh12(NucHigh12>0)); %Marker 9: nonzero Std of nucleus LowD subsection: MeanCell Image- 2*threshold (mean_cell7)< D < MeanCell Image- threshold (mean_cell6)



%--------------- Mean (Cell_Rms) - STD < Rms < Mean (Cell_Rms)
%--------------Start
%%%subplot(3,3,4);
NucHigh13=NucImage;
NucHigh13(NucHigh13>mean_cell1)=0;%% For Marker 10: Average D value of nucleus subsection: D < MeanCell Image)
%%%imshow(NucHigh2)
NucHigh13Area=NucHigh13;
NucHigh13Area(NucHigh13Area>0)=1;%For Marker 11: Area ratio of nucleus LowD subsection: D < MeanCell Image

NucHigh14=NucHigh13;
NucHigh14(NucHigh14 <= mean_cell6)=0;% For Marker 13: Average D value of nucleus subsection: MeanCell Image -LowThreshold<D < MeanCell Image)

NucHigh14Area=NucHigh14;
NucHigh14Area(NucHigh14Area>0)=1;%For Biomarker 14: Are ratio of nucleus Low D subsection: MeanCell Image - Low_Threshold < D < MeanCell Image)


HighRMSMean13= sum(sum(NucHigh13))/sum(sum(NucHigh13Area));%Marker 10: Average D value of nucleus subsection: D < MeanCell Image)
HighRMSAreaRatio13=100*(sum(sum(NucHigh13Area))/sum(sum(NucArea)));%Marker 11: Area ratio of nucleus LowD subsection: D < MeanCell Image-
HighRMSSTD13=std(NucHigh13(NucHigh13>0));% Marker 12: Nonzero Std of nucleus LowD subsection: D < MeanCell Image 

HighRMSMean14= sum(sum(NucHigh14))/sum(sum(NucHigh14Area)); %Biomarker13: Average D value of nucleus lowD subsection: MeanCell_Image - LowThreshold < D < MeanCell_Image)
HighRMSAreaRatio14=100*(sum(sum(NucHigh14Area))/sum(sum(NucArea)));%Biomarker 14: Are ratio of nucleus Low D subsection: MeanCell Image - Low_Threshold < D < MeanCell Image)

HighRMSSTD14=std(NucHigh14(NucHigh14>0)); %Biomarker 15: Nonzero Std of nucleus Low D subsection: MeanCell Image - Low_Threshold < D < MeanCell Image)



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


NucArea10Connected1= bwpropfilt(logical(NucHigh10Area),'perimeter',1);
NucArea10Connected2= bwpropfilt(logical(NucHigh10Area),'perimeter',2);
NucArea10Connected3= bwpropfilt(logical(NucHigh10Area),'perimeter',3);

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

windowSize=5;
[Marker_Edge_Ave, Marker_NearEdge_Ave, Marker_Center_Ave, Marker_Edge_Std, Marker_NearEdge_Std, Marker_Center_Std, Marker_Edge_entropy, Marker_NearEdge_entropy, Marker_Centr_entropy]= BoundryMarkers(NucImage,Mask,windowSize);
%End
%---------------------------------






% Subsection Type 1: Average D Variable in each subsection
Marker1= HighRMSMean10;%Marker 1: Average D value of nucleus subsection: D < MeanCell Image-2 threshold (mean_cell7)
Marker2=HighRMSMean11;%Marker 2: Average D value of nucleus subsection: D < MeanCell Image- threshold (mean_cell6)
Marker3=HighRMSMean13;%%Marker 3: Average D value of nucleus subsection: D < MeanCell Image (mean_cell1)
Marker4=HighRMSMean1;%%Marker 3: Average D value of nucleus subsection: D < MeanCell Image (mean_cell1)
Marker5=HighRMSMean3;%%Marker 4: Average D value of nucleus subsection: D < MeanCell Image + 1*threshold_High (mean_cell2)
Marker6=HighRMSMean5;%%Marker 5: Average D value of nucleus subsection: D < MeanCell Image + 2*threshold_High (mean_cell3)
Marker7=HighRMSMean7;%%Marker 6: Average D value of nucleus subsection: D < MeanCell Image + 3*threshold_High (mean_cell4)
Marker8=HighRMSMean9;%%Marker 6: Average D value of nucleus subsection: D < MeanCell Image + 4*threshold_High (mean_cell5)

% Subsection Type 1: Area Ratia Variable in each subsection
Marker9= HighRMSAreaRatio10;%Marker 9: Percentage Area ratio of nucleus subsection: D < MeanCell Image-2 threshold (mean_cell7)
Marker10=HighRMSAreaRatio11;%Marker 10: Percentage Area ratio of nucleus subsection: D < MeanCell Image- threshold (mean_cell6)
Marker11=HighRMSAreaRatio13;%Marker 11: Percentage Area ratio of nucleus subsection: D < MeanCell Image (mean_cell1)
Marker12=HighRMSAreaRatio1;%%Marker 12: Percentage Area ratio of nucleus subsection: D > MeanCell Image (mean_cell1)
Marker13=HighRMSAreaRatio3;%%Marker 13: APercentage Area ratio of nucleus subsection: D > MeanCell Image + 1*threshold_High (mean_cell2)
Marker14=HighRMSAreaRatio5;%%Marker 14: Percentage Area ratioof nucleus subsection: D > MeanCell Image + 2*threshold_High (mean_cell3)
Marker15=HighRMSAreaRatio7;%%Marker 15: Percentage Area ratio of nucleus subsection: D > MeanCell Image + 3*threshold_High (mean_cell4)
Marker16=HighRMSAreaRatio9;%%Marker 16: Percentage Area ratio of nucleus subsection: D > MeanCell Image + 4*threshold_High (mean_cell5)

% Subsection Type 1: Largest connected organization length (Variable) in each subsection
Marker17= 100*(sum(sum(double(NucArea10Connected1))))/sum(sum(NucArea));%Marker 17: Length of largest connected structure of nucleus subsection: D < MeanCell Image-2 threshold (mean_cell7)
Marker18= 100*(sum(sum(double(NucArea11Connected1))))/sum(sum(NucArea));%Marker 18: Length of largest connected structure of nucleus subsection: D < MeanCell Image- threshold (mean_cell6)
Marker19= 100*(sum(sum(double(NucArea13Connected1))))/sum(sum(NucArea));%Marker 19: Length of largest connected structure of nucleus subsection: D < MeanCell Image(mean_cell1)
Marker20= 100*(sum(sum(double(NucArea1Connected1))))/sum(sum(NucArea));%%Marker 20: Length of largest connected structure of nucleus subsection: D > MeanCell Image (mean_cell1)
Marker21= 100*(sum(sum(double(NucArea3Connected1))))/sum(sum(NucArea));%%Marker 21: Length of largest connected structure of nucleus subsection: D > MeanCell Image + 1*threshold_High (mean_cell2)
Marker22= 100*(sum(sum(double(NucArea5Connected1))))/sum(sum(NucArea));%%Marker 22: Length of largest connected structure of nucleus subsection: D > MeanCell Image + 2*threshold_High (mean_cell3)
Marker23= 100*(sum(sum(double(NucArea7Connected1))))/sum(sum(NucArea));%%Marker 23: Length of largest connected structure of nucleus subsection: D > MeanCell Image + 3*threshold_High (mean_cell4)
Marker24= 100*(sum(sum(double(NucArea9Connected1))))/sum(sum(NucArea));%%Marker 24: Length of largest connected structureo of nucleus subsection: D > MeanCell Image + 4*threshold_High (mean_cell5)

% Subsection Type 1: Std Variable in each subsection
Marker25= HighRMSSTD10;%Marker 25: Standard deviation of nucleus subsection: D < MeanCell Image-2 threshold (mean_cell7)
Marker26=HighRMSSTD11;%Marker 26: Standard deviation of nucleus subsection: D < MeanCell Image- threshold (mean_cell6)
Marker27=HighRMSSTD13;%Marker 27: Standard deviation of nucleus subsection: D < MeanCell Image (mean_cell1)
Marker28=HighRMSSTD1;%%Marker 28: Standard deviation of nucleus subsection: D > MeanCell Image (mean_cell1)
Marker29=HighRMSSTD3;%%Marker 29: Standard deviation of nucleus subsection: D > MeanCell Image + 1*threshold_High (mean_cell2)
Marker30=HighRMSSTD5;%%Marker 30: Standard deviation of nucleus subsection: D > MeanCell Image + 2*threshold_High (mean_cell3)
Marker31=HighRMSSTD7;%%Marker 31: Standard deviation of nucleus subsection: D > MeanCell Image + 3*threshold_High (mean_cell4)
Marker32=HighRMSSTD9;%%Marker 32: Standard deviation of nucleus subsection: D > MeanCell Image + 4*threshold_High (mean_cell5)

% Subsection Type 1: Largest connected organization length 2 (Variable) in each subsection
Marker33= 100*(sum(sum(double(NucArea10Connected2))))/sum(sum(NucArea));%Marker 8: Length of two largest connected structures of nucleus subsection: D < MeanCell Image-2 threshold (mean_cell7)
Marker34=100*(sum(sum(double(NucArea11Connected2))))/sum(sum(NucArea));%Marker 9: Length of two largest connected structures of nucleus subsection: D < MeanCell Image- threshold (mean_cell6)
Marker35=100*(sum(sum(double(NucArea13Connected2))))/sum(sum(NucArea));%Marker 9: Length of two largest connected structures of nucleus subsection: D < MeanCell Image(mean_cell1)
Marker36=100*(sum(sum(double(NucArea1Connected2))))/sum(sum(NucArea));%%Marker 10: Length of two largest connected structures of nucleus subsection: D > MeanCell Image (mean_cell1)
Marker37=100*(sum(sum(double(NucArea3Connected2))))/sum(sum(NucArea));%%Marker 11: Length of two largest connected structures of nucleus subsection: D > MeanCell Image + 1*threshold_High (mean_cell2)
Marker38=100*(sum(sum(double(NucArea5Connected2))))/sum(sum(NucArea));%%Marker 12: Length of two largest connected structures of nucleus subsection: D > MeanCell Image + 2*threshold_High (mean_cell3)
Marker39=100*(sum(sum(double(NucArea7Connected2))))/sum(sum(NucArea));%%Marker 13: Length of two largest connected structures of nucleus subsection: D > MeanCell Image + 3*threshold_High (mean_cell4)
Marker40=100*(sum(sum(double(NucArea9Connected2))))/sum(sum(NucArea));%%Marker 6: Length of two largest connected structureos of nucleus subsection: D > MeanCell Image + 4*threshold_High (mean_cell5)


% Subsection Type 1: Largest connected organization length 3 (Variable) in each subsection
Marker41= 100*(sum(sum(double(NucArea10Connected3))))/sum(sum(NucArea));%Marker 37: Length of three largest connected structures of nucleus subsection: D < MeanCell Image-2 threshold (mean_cell7)
Marker42=100*(sum(sum(double(NucArea11Connected3))))/sum(sum(NucArea));%Marker 38: Length of three largest connected structures of nucleus subsection: D < MeanCell Image- threshold (mean_cell6)
Marker43=100*(sum(sum(double(NucArea13Connected3))))/sum(sum(NucArea));%Marker 39: Length of three largest connected structures of nucleus subsection: D < MeanCell Image(mean_cell1)
Marker44=100*(sum(sum(double(NucArea1Connected3))))/sum(sum(NucArea));%%Marker 40: Length of three largest connected structures of nucleus subsection: D > MeanCell Image (mean_cell1)
Marker45=100*(sum(sum(double(NucArea3Connected3))))/sum(sum(NucArea));%%Marker 41: Length of three largest connected structures of nucleus subsection: D > MeanCell Image + 1*threshold_High (mean_cell2)
Marker46=100*(sum(sum(double(NucArea5Connected3))))/sum(sum(NucArea));%%Marker 42: Length of three largest connected structures of nucleus subsection: D > MeanCell Image + 2*threshold_High (mean_cell3)
Marker47=100*(sum(sum(double(NucArea7Connected3))))/sum(sum(NucArea));%%Marker 43: Length of three largest connected structures of nucleus subsection: D > MeanCell Image + 3*threshold_High (mean_cell4)
Marker48=100*(sum(sum(double(NucArea9Connected3))))/sum(sum(NucArea));%%Marker 44: Length of three largest connected structureos of nucleus subsection: D < MeanCell Image + 4*threshold_High (mean_cell5)

Marker49=TraditionalRms;%Traditional AveRms of Nucleus


%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
% Subsection Type 2: Average D Variable in each subsection
Marker50= HighRMSMean12;%Marker 46: Average D value of nucleus subsection: MeanCell Image-2 threshold (mean_cell7)<D < MeanCell Image-1 threshold (mean_cell6)
Marker51=HighRMSMean14;%Marker 47: Average D value of nucleus subsection: MeanCell Image- threshold (mean_cell6)<D < MeanCell Image (mean_cell1)
Marker52=HighRMSMean2;%%Marker 48: Average D value of nucleus subsection: MeanCell Image (mean_cell1)<D < MeanCell Image + 1*threshold_High (mean_cell2))
Marker53=HighRMSMean4;%%Marker 49: Average D value of nucleus subsection: MeanCell Image + 1*threshold_High (mean_cell2)<D < MeanCell Image + 2*threshold_High (mean_cell3)
Marker54=HighRMSMean6;%%Marker 50: Average D value of nucleus subsection: MeanCell Image + 2*threshold_High (mean_cell3)<D < MeanCell Image + 3*threshold_High (mean_cell4)
Marker55=HighRMSMean8;%%Marker 51: Average D value of nucleus subsection: MeanCell Image + 3*threshold_High (mean_cell4)<D < MeanCell Image + 3*threshold_High (mean_cell5)
Marker56=HighRMSMean9;%%Marker 52: Average D value of nucleus subsection: MeanCell Image + 4*threshold_High (mean_cell5)<D  


% Subsection Type 2: AriaRatio Variable in each subsection
Marker57= HighRMSAreaRatio12;%Marker 53: Percentage Area ratio of nucleus subsection: MeanCell Image-2 threshold (mean_cell7)<D < MeanCell Image-1 threshold (mean_cell6)
Marker58=HighRMSAreaRatio14;%Marker 54: Percentage Area ratio of nucleus subsection: MeanCell Image- threshold (mean_cell6)<D < MeanCell Image (mean_cell1)
Marker59=HighRMSAreaRatio2;%%Marker 55: Percentage Area ratio of nucleus subsection: MeanCell Image (mean_cell1)<D < MeanCell Image + 1*threshold_High (mean_cell2))
Marker60=HighRMSAreaRatio4;%%Marker 56: Percentage Area ratio of nucleus subsection: MeanCell Image + 1*threshold_High (mean_cell2)<D < MeanCell Image + 2*threshold_High (mean_cell3)
Marker61=HighRMSAreaRatio6;%%Marker 57: Percentage Area ratio of nucleus subsection: MeanCell Image + 2*threshold_High (mean_cell3)<D < MeanCell Image + 3*threshold_High (mean_cell4)
Marker62=HighRMSAreaRatio8;%%Marker 58: Percentage Area ratio of nucleus subsection: MeanCell Image + 3*threshold_High (mean_cell4)<D < MeanCell Image + 3*threshold_High (mean_cell5)
Marker63=HighRMSAreaRatio9;%%Marker 59: Percentage Area ratio of nucleus subsection: MeanCell Image + 4*threshold_High (mean_cell5)<D  

% Subsection Type 2: Largest connected organization length 2 (Variable) in each subsection
Marker64= 100*(sum(sum(double(NucArea12Connected1))))/sum(sum(NucArea));%Marker 58: Length of largest connected structure of nucleus subsection: MeanCell Image-2 threshold (mean_cell7)<D < MeanCell Image-1 threshold (mean_cell6)
Marker65=100*(sum(sum(double(NucArea14Connected1))))/sum(sum(NucArea));%Marker 59: Length of largest connected structure of nucleus subsection: MeanCell Image- threshold (mean_cell6)<D < MeanCell Image (mean_cell1)
Marker66=100*(sum(sum(double(NucArea2Connected1))))/sum(sum(NucArea));%%Marker 60: Length of largest connected structure of nucleus subsection: MeanCell Image (mean_cell1)<D < MeanCell Image + 1*threshold_High (mean_cell2))
Marker67=100*(sum(sum(double(NucArea4Connected1))))/sum(sum(NucArea));%%Marker 61: Length of largest connected structure of nucleus subsection: MeanCell Image + 1*threshold_High (mean_cell2)<D < MeanCell Image + 2*threshold_High (mean_cell3)
Marker68=100*(sum(sum(double(NucArea6Connected1))))/sum(sum(NucArea));%%Marker 62: Length of largest connected structure of nucleus subsection: MeanCell Image + 2*threshold_High (mean_cell3)<D < MeanCell Image + 3*threshold_High (mean_cell4)
Marker69=100*(sum(sum(double(NucArea8Connected1))))/sum(sum(NucArea));%%Marker 63: Length of largest connected structure of nucleus subsection: MeanCell Image + 3*threshold_High (mean_cell4)<D < MeanCell Image + 3*threshold_High (mean_cell5)
Marker70=100*(sum(sum(double(NucArea9Connected1))))/sum(sum(NucArea));%%Marker 64: Length of largest connected structure of nucleus subsection: MeanCell Image + 4*threshold_High (mean_cell5)<D  


% Subsection Type 2: Standard deviation Variable in each subsection
Marker71= HighRMSSTD12;%Marker 65: Standard deviation of nucleus subsection: MeanCell Image-2 threshold (mean_cell7)<D < MeanCell Image-1 threshold (mean_cell6)
Marker72=HighRMSSTD14;%Marker 66: Standard deviation of nucleus subsection: MeanCell Image- threshold (mean_cell6)<D < MeanCell Image (mean_cell1)
Marker73=HighRMSSTD2;%%Marker 67: Standard deviation of nucleus subsection: MeanCell Image (mean_cell1)<D < MeanCell Image + 1*threshold_High (mean_cell2))
Marker74=HighRMSSTD4;%%Marker 68: Standard deviation of nucleus subsection: MeanCell Image + 1*threshold_High (mean_cell2)<D < MeanCell Image + 2*threshold_High (mean_cell3)
Marker75=HighRMSSTD6;%%Marker 69: Standard deviation of nucleus subsection: MeanCell Image + 2*threshold_High (mean_cell3)<D < MeanCell Image + 3*threshold_High (mean_cell4)
Marker76=HighRMSSTD8;%%Marker 70: Standard deviation of nucleus subsection: MeanCell Image + 3*threshold_High (mean_cell4)<D < MeanCell Image + 3*threshold_High (mean_cell5)
Marker77=HighRMSSTD9;%%Marker 71: Standard deviation of nucleus subsection: MeanCell Image + 4*threshold_High (mean_cell5)<D  

% Subsection Type 2: two Largest connected organization length 2 (Variable) in each subsection
Marker78= 100*(sum(sum(double(NucArea12Connected2))))/sum(sum(NucArea));%Marker 72: Length of two largest connected structure of nucleus subsection: MeanCell Image-2 threshold (mean_cell7)<D < MeanCell Image-1 threshold (mean_cell6)
Marker79=100*(sum(sum(double(NucArea14Connected2))))/sum(sum(NucArea));%Marker 73: Length of two largest connected structure of nucleus subsection: MeanCell Image- threshold (mean_cell6)<D < MeanCell Image (mean_cell1)
Marker80=100*(sum(sum(double(NucArea2Connected2))))/sum(sum(NucArea));%%Marker 74: Length of two largest connected structure of nucleus subsection: MeanCell Image (mean_cell1)<D < MeanCell Image + 1*threshold_High (mean_cell2))
Marker81=100*(sum(sum(double(NucArea4Connected2))))/sum(sum(NucArea));%%Marker 75: Length of two largest connected structure of nucleus subsection: MeanCell Image + 1*threshold_High (mean_cell2)<D < MeanCell Image + 2*threshold_High (mean_cell3)
Marker82=100*(sum(sum(double(NucArea6Connected2))))/sum(sum(NucArea));%%Marker 76: Length of  twolargest connected structure of nucleus subsection: MeanCell Image + 2*threshold_High (mean_cell3)<D < MeanCell Image + 3*threshold_High (mean_cell4)
Marker83=100*(sum(sum(double(NucArea8Connected2))))/sum(sum(NucArea));%%Marker 77: Length of two largest connected structure of nucleus subsection: MeanCell Image + 3*threshold_High (mean_cell4)<D < MeanCell Image + 3*threshold_High (mean_cell5)
Marker84=100*(sum(sum(double(NucArea9Connected2))))/sum(sum(NucArea));%%Marker 78: Length of two largest connected structure of nucleus subsection: MeanCell Image + 4*threshold_High (mean_cell5)<D  

% Subsection Type 2: two Largest connected organization length 2 (Variable) in each subsection
Marker85= 100*(sum(sum(double(NucArea12Connected3))))/sum(sum(NucArea));%Marker 79: Length of two largest connected structure of nucleus subsection: MeanCell Image-2 threshold (mean_cell7)<D < MeanCell Image-1 threshold (mean_cell6)
Marker86=100*(sum(sum(double(NucArea14Connected3))))/sum(sum(NucArea));%Marker 80: Length of two largest connected structure of nucleus subsection: MeanCell Image- threshold (mean_cell6)<D < MeanCell Image (mean_cell1)
Marker87=100*(sum(sum(double(NucArea2Connected3))))/sum(sum(NucArea));%%Marker 81: Length of two largest connected structure of nucleus subsection: MeanCell Image (mean_cell1)<D < MeanCell Image + 1*threshold_High (mean_cell2))
Marker88=100*(sum(sum(double(NucArea4Connected3))))/sum(sum(NucArea));%%Marker 82: Length of two largest connected structure of nucleus subsection: MeanCell Image + 1*threshold_High (mean_cell2)<D < MeanCell Image + 2*threshold_High (mean_cell3)
Marker89=100*(sum(sum(double(NucArea6Connected3))))/sum(sum(NucArea));%%Marker 83: Length of  twolargest connected structure of nucleus subsection: MeanCell Image + 2*threshold_High (mean_cell3)<D < MeanCell Image + 3*threshold_High (mean_cell4)
Marker90=100*(sum(sum(double(NucArea8Connected3))))/sum(sum(NucArea));%%Marker 84: Length of two largest connected structure of nucleus subsection: MeanCell Image + 3*threshold_High (mean_cell4)<D < MeanCell Image + 3*threshold_High (mean_cell5)
Marker91=100*(sum(sum(double(NucArea9Connected3))))/sum(sum(NucArea));%%Marker 85: Length of two largest connected structure of nucleus subsection: MeanCell Image + 4*threshold_High (mean_cell5)<D  


%Differential Markers
%LowD subsections Differential Markers(D)
Marker92=HighRMSMean1-HighRMSMean13;
Marker93=HighRMSMean13-HighRMSMean11;
Marker94=HighRMSMean11-HighRMSMean10;
%HighD subsections Differential Markers(D)

Marker95=HighRMSMean3-HighRMSMean1;
Marker96=HighRMSMean5-HighRMSMean3;
Marker97=HighRMSMean7-HighRMSMean5;
Marker98=HighRMSMean9-HighRMSMean7;

%LowD subsections Differential Markers(AreaRatio)
Marker99=abs(HighRMSAreaRatio10-HighRMSAreaRatio11);
Marker100=abs(HighRMSAreaRatio11-HighRMSAreaRatio13);
Marker101=abs(HighRMSAreaRatio1-HighRMSAreaRatio13);
%HighD subsections Differential Markers(AreaRatio)
Marker102=abs(HighRMSAreaRatio1-HighRMSAreaRatio3);
Marker103=abs(HighRMSAreaRatio3-HighRMSAreaRatio5);
Marker104=abs(HighRMSAreaRatio5-HighRMSAreaRatio7);
Marker105=abs(HighRMSAreaRatio7-HighRMSAreaRatio9);

%HighD subsections



Marker106=std([Marker92 Marker93 Marker94]); %Std of differential Rms-markers in LowD subsection
Marker107=std([Marker95 Marker96 Marker97 Marker98]);%Std of differential Rms-markers in HighD subsection
Marker108=std([Marker92 Marker93 Marker94 Marker95 Marker96 Marker97 Marker98]);%Std of differential Rms-markers across all subsection (Low and High D)

Marker109=std([Marker99 Marker100 Marker101]);%Std of differential AreaRatio-markers in LowD subsection
Marker110=std([Marker102 Marker103 Marker104 Marker105]);%Std of differential AreaRatio-markers in HighD subsection
Marker111=std([Marker99 Marker100 Marker101 Marker102 Marker103 Marker104 Marker105]);%Std of differential AreaRatio-markers across all subsection (Low and High D)




%%%%Boundry Markers
Marker112=Marker_Edge_Ave;% Edge Ave Rms
Marker113=Marker_NearEdge_Ave;% Near Edge Ave Rms
Marker114=Marker_Center_Ave;% Center  Ave Rms

Marker115=Marker_Edge_Std;
Marker116=Marker_NearEdge_Std;
Marker117=Marker_Center_Std;

Marker118=Marker_Edge_entropy;
Marker119=Marker_NearEdge_entropy;
Marker120=Marker_Centr_entropy;

Marker121=Marker_Edge_Ave/Marker_NearEdge_Ave;
Marker122=Marker_Edge_Ave/Marker_Center_Ave;
Marker123=Marker_NearEdge_Ave/Marker_Center_Ave;

Marker124=Marker_Edge_Ave*Marker_NearEdge_Ave;
Marker125=Marker_Edge_Ave*Marker_Center_Ave;
Marker126=Marker_NearEdge_Ave*Marker_Center_Ave;



MakerArray=[Marker1 Marker2 Marker3 Marker4 Marker5 Marker6 Marker7 Marker8 Marker9 Marker10 Marker11 Marker12 Marker13 Marker14 Marker15 Marker16 Marker17 Marker18 Marker19 Marker20 Marker21 Marker22 Marker23 Marker24 Marker25 Marker26 Marker27 Marker28 Marker29 Marker30 Marker31 Marker32 Marker33 Marker34 Marker35 Marker36 Marker37 Marker38 Marker39 Marker40 Marker41 Marker42 Marker43 Marker44 Marker45 Marker46 Marker47 Marker48 Marker49 Marker50 Marker51 Marker52 Marker53 Marker54 Marker55 Marker56 Marker57 Marker58 Marker59 Marker60 Marker61 Marker62 Marker63 Marker64 Marker65 Marker66 Marker67 Marker68 Marker69 Marker70 Marker71 Marker72 Marker73 Marker74 Marker75 Marker76 Marker77 Marker78 Marker79 Marker80 Marker81 Marker82 Marker83 Marker84 Marker85 Marker86 Marker87 Marker88 Marker89 Marker90 Marker91 Marker92 Marker93 Marker94 Marker95 Marker96 Marker97 Marker98 Marker99 Marker100 Marker101 Marker102 Marker103 Marker104 Marker105 Marker106 Marker107 Marker108 Marker109 Marker110 Marker111 Marker112 Marker113 Marker114 Marker115 Marker116 Marker117 Marker118 Marker119 Marker120 Marker121 Marker122 Marker123 Marker124 Marker125 Marker126];

%Questions can be adressed to alid@northwestern.edu
