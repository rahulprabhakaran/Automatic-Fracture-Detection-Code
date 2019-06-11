clc
clear all
close all
format compact

%% select files to be opened
    disp('Select multiple Ridges files (they must be in a single folder)');
    [InFileListShort, pathname] = uigetfile('*.jpg;*.tif;*.png','Select input image files','MultiSelect','on');
    
    if not(iscell(InFileListShort))
       InFileListShort = {InFileListShort};
    end
    
    replicatePath = repmat(cellstr(pathname),size(InFileListShort));
    InFileList = strcat(replicatePath,InFileListShort);    
           
    clear InFileListShort pathname replicatePath;
    funPath = fileparts(which('Segmentation_Otsu_Thresholding.m'));
    addpath(genpath(funPath));   
    
    Path_Segmented_Ridges = 'D:\PhD\Automatic_Detection\Probabilistic_Edges\Ortho_1_107\Ortho_107_New\Segmented_Ridges';
    Path_Segmented_Ridges_Overlay = 'D:\PhD\Automatic_Detection\Probabilistic_Edges\Ortho_1_107\Ortho_107_New\Segmented_Ridges_Overlay';
    Path_Skeletons = 'D:\PhD\Automatic_Detection\Probabilistic_Edges\Ortho_1_107\Ortho_107_New\Skeletons';
    Path_Fitted_Curves = 'D:\PhD\Automatic_Detection\Probabilistic_Edges\Ortho_1_107\Ortho_107_New\Fitted_Curves';

%%
coordRefSysCode = 32632;
for m = 1:size(InFileList, 2)
    tic   
    % counter
    disp(' ');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(['This is image ' num2str(m) ' out of ' num2str(size(InFileList, 2))]);
    disp(['Performing Otsu Thresholding Segmentation for Image  ' num2str(m) ' out of ' num2str(size(InFileList, 2))]);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    
    % Ridges image files
    currentImageFileName = InFileList{m};
    [currentPath,currentImageShortName,~] = fileparts(currentImageFileName);
    
    % Orthomosaics image files
    % change the number of characters to be subtracted so that the 
    % currentImageShortName2 character array corresponds to the name of the
    % original image files
    currentPath2 = [currentPath(1:length(currentPath)) '\'];
%     currentImageShortName2 = [currentImageShortName(1:length(currentImageShortName)) '.tif'];
    currentImageShortName2 = [currentImageShortName(8:length(currentImageShortName)) '.tif'];
    
    % open Ridges image and the original tiled orthomosaic
    imageIN = imread(currentImageFileName); % imageIN is read
    %imageIN = imresize(imageIN,[1001 1021]);
%     imageIN = imageIN(:,:,1);
%     [~,R] = geotiffread(InFileList{i});
    
%     imageIN = imrotate(imageIN,180);
%     imageIN_1 = rgb2gray(imread([currentPath2 currentImageShortName2]));       % rgb version is used for overlay
      imageIN_1 = imread([currentPath2 currentImageShortName2]);                 % if image is grayscale use this statement
      %imageIN_1 = imresize(imageIN_1,[1000 1000]); 
      %imageIN_1 = rgb2gray(imageIN_1); 
      imageIN_4 = imageIN;
%     imageIN_1 = imrotate(imageIN_1,180);
    
    % resizing original image for overlay
%     imageIN_1 = imresize(imageIN_1,0.25);
    
    % estimating background using morphological opening
    background = imopen(imageIN,strel('disk',5));
    %imshow(background)
    
    % subtracting the background image from the original image
    imageIN_2 = imageIN - background;
    %imshow(imageIN_2)
    
    % increasing the image contrast
    imageIN_3 = imadjust(imageIN_2);
    %imshow(imageIN_3)

    % thresholding the image, set number of pixels to be removed, P = 10
    bw = imbinarize(imageIN_3);
    bw = bwareaopen(bw, 20);
    %imshow(bw)
    
    % identifying objects within the image
    cc = bwconncomp(bw, 4);
    
%     % examining one object
%     grain = false(size(bw));
%     grain(cc.PixelIdxList{1049}) = true;
    
    % creating labels
    labeled = labelmatrix(cc);
    whos labeled;
    RGB_label = label2rgb(labeled, @spring, 'c', 'shuffle');
    
    % removing isolated clusters
    Pixel_List = cc.PixelIdxList';
    for j=1:length(Pixel_List)
     Pixel_List_length(j,1)=length(Pixel_List{j,1});
    end
    imageIN_4 = labeled;
    %imshow(imageIN_4)

    Isolation_Threshold = 10;
    for j=1:length(Pixel_List)
     if Pixel_List_length(j,1)< Isolation_Threshold
      imageIN_4(imageIN_4==j)=0;
     end
    end
    %imshow(imageIN_4)
    
    j=1;
    tic
    for k=1:length(Pixel_List)
     [r,c]=find(imageIN_4==k);
     if isempty(r)~=1  
      Cluster_List{j,1}=[r c];
      j=j+1;
     end
    end
    toc
    
%     set(gca,'Ydir','normal')
%     for i=1:length(Cluster_List)
%        set(gca,'Ydir','reverse')
%        scatter(Cluster_List{i,1}(:,1),Cluster_List{i,1}(:,2),'Filled','b');
%        hold on
%     end    
    
   %  converting to image
   imageIN_4(imageIN_4>0)=255;
%    imageIN_4=mat2gray(imageIN_4);
  imageIN_1=double(255 * mat2gray(imageIN_1));
  imageIN_4=double(255 * mat2gray(imageIN_4));
  %imshow(imageIN_1)

   % 
%     imageIN_4=im2uint8(imageIN_4);
   
%  writing the segmented ridges and its overlay to the output folders
%    outSuffix_Segmented_Ridges = [currentImageShortName(1:length(currentImageShortName)-7) '_Segmented_Ridges' '.tif'];
%    outSuffix_Segmented_Ridges_Overlay = [currentImageShortName(1:length(currentImageShortName)-7) '_Segmented_Ridges_Overlay' '.tif'];
   
   outSuffix_Segmented_Ridges = [currentImageShortName '_Segmented_Ridges' '.tif'];
   outSuffix_Segmented_Ridges_Overlay = [currentImageShortName '_Segmented_Ridges_Overlay' '.tif'];

   OutFileName = [Path_Segmented_Ridges '\' outSuffix_Segmented_Ridges];
   disp(['Writing the Segmented Ridges for Image ' num2str(m) ' out of ' num2str(size(InFileList, 2))]);
  imwrite(imageIN_4, OutFileName); 
%    geotiffwrite(OutFileName,imageIN_4,R,'CoordRefSysCode', coordRefSysCode);

   %imageIN_1 = imresize(imageIN_1 ,[1000 1000]);
   overlay_plus_segmented_ridges=CSHRMgetOverlay(imageIN_1,imageIN_4);
   OutFileName = [Path_Segmented_Ridges_Overlay '\' outSuffix_Segmented_Ridges_Overlay];
   disp(['Writing the Segmented Ridges and Overlay for Image ' num2str(m) ' out of ' num2str(size(InFileList, 2))]);
   imwrite(overlay_plus_segmented_ridges, OutFileName); 
%   geotiffwrite(OutFileName,overlay_plus_segmented_ridges,R,'CoordRefSysCode', coordRefSysCode);
   
   
   % Skeletonizing the segmented ridges using Otsu Thresholding and
   % calculating branch points and end points for each cluster
   skelImg   = bwmorph(imageIN_4, 'thin', 'inf');
   branchImg = bwmorph(skelImg, 'branchpoints');
   endImg    = bwmorph(skelImg, 'endpoints');

   [row, column] = find(endImg);
   endPts        = [row column];

   [row, column] = find(branchImg);
   branchPts     = [row column];
   cNumBranchPoints = length(branchPts);
   
   %  writing the skeletonized ridges to the output folders
%    outSuffix_Skeletons = [currentImageShortName(1:length(currentImageShortName)-7) '_Skeletons' '.tif'];
   outSuffix_Skeletons = [currentImageShortName '_Skeletons' '.tif'];
   OutFileName = [Path_Skeletons '\' outSuffix_Skeletons];
   disp(['Writing the Skeletons for Image ' num2str(m) ' out of ' num2str(size(InFileList, 2))]);
   
   %skelImg2 = imresize(skelImg,4);
   skelImg2 = mat2gray(skelImg);
   
   imwrite(~skelImg2, OutFileName); 
   %geotiffwrite(OutFileName,~skelImg2,R,'CoordRefSysCode', coordRefSysCode);
   
   
   % concatenating the clusters into one large matrix
    Clusters=Cluster_List{1,1};
    for i=2:length(Cluster_List)
     A= Cluster_List{i,1};
     Clusters=[Clusters;A];
    end

   % finding the branch points and end points associated with each cluster
   % and storing them in two additional columns
   % clusters which have no branches and end points are identified and
   % stored in no_branches and no_ends arrays
   
   y=1;
   z=1;
   for i=1:length(Cluster_List)
    Cluster_List{i,2}=intersect(branchPts,Cluster_List{i,1},'rows') ;
    Cluster_List{i,3}=intersect(endPts,Cluster_List{i,1},'rows') ;
     if isempty(Cluster_List{i,2})==1
      no_branches(y,1) = i;
      y=y+1;
     end
     if isempty(Cluster_List{i,3})==1
      no_ends(z,1) = i;
      z=z+1;
     end   
   end

   % removing clusters that have no end points. It seems these clusters
   % are very close to existing clusters and escape the bwmorph function
   if exist('Nodes','var')==1
    length_Cluster_List = length(Cluster_List);
    idx = find(no_ends==length_Cluster_List); 
    Cluster_List(no_ends(1:idx),:)=[];
   end 
   
   % fitting curves through the clusters using function from David Legland.
   % compute coeffs of each individual branch and returns a matrix of
   % labels for each fitted curve. I use this function for the time being.
   % It does not use the endPts and branchPts calculated using the
   % bwmorph call but calculates end points and branch points using bwlabel,
   % bwconncomp and regionprops functions.
   skelImg = imrotate(skelImg,180);
   [coeffs, curve_matrix] = polynomialCurveSetFit(skelImg, 5);
   
     
   %writing the polynomial fit to a table and stored in the output folder
   if isempty(coeffs)~=1
     
    % Obtaining the polynomial points for each fitted curve using coeffs for
     %each curve
    for n = 1:length(coeffs)
     Poly_Points{n,1}= drawPolynomialCurve([0 1], coeffs{n});
    end
     
     
    %figure(9); imshow(~skelImg); hold on;
%         for i = 1:length(coeffs)
%            Poly_Points{i,1}= drawPolynomialCurve([0 1], coeffs{i});
%    
%             set(hc, 'linewidth', 2, 'color', 'r');
%         end
  
          for i=1:length(Poly_Points)
          Poly_Points_Table_Header {i}= (['Polyline_' ,num2str(i)]);
        end
    Poly_Points_Table=cell2table(Poly_Points);
       
    outSuffix_Fitted_Curves = [currentImageShortName(1:length(currentImageShortName)-7) '_Fitted_Curves' '.mat'];
    outSuffix_Fitted_Curves = [currentImageShortName '_Fitted_Curves' '.mat'];
    OutFileName = [Path_Fitted_Curves '\' outSuffix_Fitted_Curves];
    disp(['Writing Polyline Points for Image ' num2str(m) ' out of ' num2str(size(InFileList, 2))]);
    save(OutFileName,'Poly_Points_Table');
   
   else
      disp(['Image ' num2str(i) ' is empty. No Polylines to write' ]); 
   writetable(Poly_Points_Table,OutFileName);
   end
   
   clearvars coeffs; 
   clearvars Poly_Points; 
   clearvars Poly_Points_Table; 
   clearvars curve_matrix; 
   clearvars n; 
   clearvars R
   clearvars Cluster_List
   clearvars no_ends
   
   toc
end

%%
% the code below was written to test the efficiency of Otsu thresholding
% using a single image. The loop above performs the same for 

%% reading the image
% I = imread('D:\PhD\Automatic_Detection\Tiled_Orthomosaics\ortho_1\Ridges\Ortho_1_107_Ridges.tif');
% figure(1)
% imshow(I)
% 
% %% Using morphological opening to estimate the background
% 
% background = imopen(I,strel('disk',15));
% 
% % Display the Background Approximation as a Surface
% figure
% surf(double(background(1:8:end,1:8:end))),zlim([0 255]);
% ax = gca;
% ax.YDir = 'reverse';
% 
% %% subtracting the background image from the original image
% I2 = I - background;
% figure(2)
% imshow(I2)
% 
% %% increasing the image contrast
% figure(3)
% I3 = imadjust(I2);
% imshow(I3);
% 
% %% thresholding the image
% bw = imbinarize(I3);
% bw = bwareaopen(bw, 50);
% figure(4)
% imshow(bw)
% 
% %% identifying objects within the image
% cc = bwconncomp(bw, 4);
% 
% %% examining one object
% grain = false(size(bw));
% grain(cc.PixelIdxList{1049}) = true;
% imshow(grain);
% 
% %% view all objects
% labeled = labelmatrix(cc);
% whos labeled
% RGB_label = label2rgb(labeled, @spring, 'c', 'shuffle');
% imshow(RGB_label)
% 
% %% removing isolated clusters
% Pixel_List = cc.PixelIdxList';
% for i=1:length(Pixel_List)
%     Pixel_List_length(i,1)=length(Pixel_List{i,1});
% end
% I4 = labeled;
% 
% 
% Isolation_Threshold = 10;
% for i=1:length(Pixel_List)
%    if Pixel_List_length(i,1)< Isolation_Threshold
%       I4(I4==i)=0;
%    end
% end
% 
% %%
% I4(I4>0)=255;
% I4=mat2gray(I4);
% imshow(I4);
% 
% %%
% I5 = imread('D:\PhD\Automatic_Detection\Tiled_Orthomosaics\ortho_1\Ortho_1_107.png');
% I5=rgb2gray(I5);
% imshow(I5);
% %%
% overlay_plus_segmented_ridges=CSHRMgetOverlay(I5,I4);
% imshow(overlay_plus_segmented_ridges);
% 
% %%
% imwrite(I4, 'D:\PhD\Automatic_Detection\Tiled_Orthomosaics\ortho_1\Segmented_Ridges\Ortho_1_107_Segmented_Ridges.png')

