%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This script converts ridge ensemble images into binary images based on
%  a ridge threshold. The ridge ensemble image is first normalized from
%  grayscale (0-255) to a range from (0-1). A non-linear sigmoid function
%  is then applied on the normalized image, so that picking the ridge 
%  intensity threshold is easier. Once the threshold is set, then all ridge
%  pixels above the threshold is set to 1 and the remaining pixels are set
%  to zero. The user has to compare the thresholded ridges with the original 
%  image for a confident assessment of high probability ridges. If certain
%  fractures (ridges) are missed in the thresholded ridge image, then the
%  shearlet parameter selection needs to be revisited. If the thresholded
%  ridge image is satisfactory, then it can be saved as binarized image
%
%  Summary of Operations  
%   ridge ensemble -> normalized ridge ensemble -> sigmoided normalized...
%  ..ridge ensemble -> binarized ridge

clc
clear all
close all
format compact

%% select ridge ensemble images to convert to binary imaages

    disp('Select multiple image files (they must be in a single folder)');
    [InFileListShort, pathname] = uigetfile('*.jpg;*.tif;*.png','Select input image files','MultiSelect','on');
    
    if not(iscell(InFileListShort))
       InFileListShort = {InFileListShort};
    end
    
    replicatePath = repmat(cellstr(pathname),size(InFileListShort));
    InFileList = strcat(replicatePath,InFileListShort);    
           
    clear pathname replicatePath;
    funPath = fileparts(which('W125_G63_SPO2_SL3_AL0.5_OCT_3.5_MC00.png'));
    addpath(genpath(funPath));  
    
%%
% reading the image
C_Ridges = imread(InFileList{1});

% displaying the ridge ensemble
figure(1)
imshow(C_Ridges)

% displaying the normalized image
C_Ridges_norm = mat2gray(C_Ridges);
figure(2)
h=imshow(C_Ridges_norm);

% using sigmoid nonlinearity on the normalized image and diaplaying it
[m,n]=size(C_Ridges_norm);
for i=1:m
  for j=1:n  
    if   C_Ridges_norm(i,j)~=0
      C_Ridges_norm_sigmoid (i,j) = 1 / (1 + exp((-1)*C_Ridges_norm(i,j)));
    end  
  end
  i
end
figure(3)
imshow(C_Ridges_norm_sigmoid)

%% Pick the value of threshold and observe the change in ridges

C_Ridges_norm_thresh = C_Ridges_norm_sigmoid;
threshold=0.51;
C_Ridges_norm_thresh(C_Ridges_norm_thresh<threshold)=0;
C_Ridges_norm_thresh(C_Ridges_norm_thresh>threshold)=1;
C_Ridges_norm_thresh= im2uint8(C_Ridges_norm_thresh);
figure(4)
imshow(C_Ridges_norm_thresh)
title(['Normalized Ridge Intensity Threshold: ',num2str(threshold)])

%% after the threshold is decided, run the following segment to create
%  binarized image and save it

% specify where file has to be saved
outfolder ='D:\Github_Test\P_Ridges\';
imwrite(C_Ridges_norm_thresh,strcat(outfolder,'Bin_',InFileListShort{1}));

%% if there are multiple images, for which a single threshold is to be applied

for i=1:length(InFileList)
 imageIN = imread(InFileList{i});
 [imageIN_Xpixels(i,1),imageIN_Ypixels(i,1)]=size(imageIN);
 clearvars imageIN
end

output_folder ='D:\Github_Test\P_Ridges\';

for k=1:length(InFileList)
   tic
   m=imageIN_Xpixels(k,1);
   n=imageIN_Ypixels(k,1);
   C_Ridges = imread(InFileList{k}) ;
   C_Ridges_norm = mat2gray(C_Ridges);
   for i=1:m
    for j=1:n
        if   C_Ridges_norm(i,j)~=0
          C_Ridges_norm_sigmoid (i,j) = 1 / (1 + exp((-1)*C_Ridges_norm(i,j)));
        end  
    end
   end   
    threshold=0.55;
    C_Ridges_norm_sigmoid_thresh = C_Ridges_norm_sigmoid;
    C_Ridges_norm_sigmoid_thresh(C_Ridges_norm_sigmoid_thresh<threshold)=0;
    C_Ridges_norm_sigmoid_thresh(C_Ridges_norm_sigmoid_thresh>threshold)=1;
    imwrite(C_Ridges_norm_sigmoid_thresh,strcat(output_folder,'Bin_',InFileListShort{k}));
    disp(['Writing the binarized ridge ensemble for image: ',num2str(k)])
    toc
    clearvars m n C_Ridges C_Ridges_norm C_Ridges_norm_sigmoid C_Edges_norm_sigmoid_thresh
end    
