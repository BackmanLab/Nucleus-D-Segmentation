function [Marker_Edge_Ave, Marker_NearEdge_Ave, Marker_Center_Ave, Marker_Edge_Std, Marker_NearEdge_Std, Marker_Center_Std, Marker_Edge_entropy, Marker_NearEdge_entropy, Marker_Center_entropy]= BoundryMarkers(NucImage,Mask,windowSize)
clear Boundry blurredBoundry kernel CenterBoundry B BB blurredBoundry1 blurredBoundry2 Center
B = bwboundaries(Mask,'noholes') ;
BB=cell2mat(B);
[m_nuc,n_nuc]=size(NucImage);
Boundry(m_nuc,n_nuc)=0;
blurredBoundry1(m_nuc,n_nuc)=0;
blurredBoundry2(m_nuc,n_nuc)=0;
CenterBoundry(m_nuc,n_nuc)=0;
Center(m_nuc,n_nuc)=0;
[m, n]=size(BB);

for i=1:m
 Boundry(BB(i,1),BB(i,2))=1;
end


%windowSize = 6; % Whatever you want.
kernel = ones(windowSize, windowSize) / windowSize ^ 2;
blurredBoundry1 = imfilter(Boundry, kernel, 'symmetric');
blurredBoundry1(blurredBoundry1>0)=1;
blurredBoundry1(NucImage==0)=0;


%Correction_Coef=(sum(sum(blurredBoundry1))+sum(sum(Boundry)))/(2*sum(sum(blurredBoundry1)));
Image_Bundry1=NucImage;
Image_Bundry1(blurredBoundry1==0)=0;
Image_Bundry1(Mask==0)=0;

%------Marker Second circle

windowSize = windowSize*2; % Whatever you want.
kernel = ones(windowSize, windowSize) / windowSize ^ 2;
blurredBoundry2 = imfilter(Boundry, kernel, 'symmetric');
blurredBoundry2(blurredBoundry2>0)=1;
blurredBoundry2(NucImage==0)=0;

Image_Bundry2=NucImage;
Image_Bundry2(blurredBoundry2==0)=0;
Image_Bundry2=Image_Bundry2-Image_Bundry1;
Image_Bundry2(Mask==0)=0;

%%---Marker Center
Max_MaskHight=max(max(BB(:,1)))-min(min(BB(:,1)));
Max_Maskwidth=max(max(BB(:,2)))-min(min(BB(:,2)));


c1=double(round(mean(BB(:,1))));
c2=double(round(mean(BB(:,2))));
Center(c1(1),c2(1))=1;
windowSize = round(0.1*min(abs([Max_MaskHight Max_Maskwidth])),0); % The size of disk is correlated with size of nucleus.
kernel = fspecial('disk',windowSize);
kernel(kernel>0)=1;
CenterBoundry = imfilter(Center, kernel, 'symmetric');
CenterBoundry(CenterBoundry>0)=1;

Image_Center=NucImage;
Image_Center(CenterBoundry==0)=0;

Marker_Edge_Ave=(sum(sum(Image_Bundry1))/sum(sum(blurredBoundry1)));
Marker_NearEdge_Ave=(sum(sum(Image_Bundry2))/sum(sum(blurredBoundry2)));
Marker_Center_Ave=sum(sum(Image_Center))/sum(sum(CenterBoundry));

Marker_Edge_Std=std(Image_Bundry1(Image_Bundry1>0));
Marker_NearEdge_Std=std(Image_Bundry2(Image_Bundry2>0));
Marker_Center_Std=std(Image_Center(Image_Center>0));

Marker_Edge_entropy = entropy(double(Image_Bundry1));
Marker_NearEdge_entropy = entropy(double(Image_Bundry2));
Marker_Center_entropy = entropy(double(Image_Center));


%imshow(CenterBoundry*255)

%rgbImage1 = cat(3, 2*Image_Bundry1, 0.5*Image_Bundry1, 0.5*Image_Bundry1);
%rgbImage2 = cat(3, 0.5*Image_Bundry2, 1*Image_Bundry2, 2*Image_Bundry2);
%rgbImage3 = cat(3, 1*Image_Center, 2*Image_Center, 0.5*Image_Center);
%RGBIMAGE=2*rgbImage1+2*rgbImage2+2*rgbImage3;
%imshow(RGBIMAGE)
