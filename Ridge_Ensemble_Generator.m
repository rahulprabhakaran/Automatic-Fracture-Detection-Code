%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Automatic Fracture Trace Detection with Complex Shearlets
%  Copyright Rahul Prabhakaran, TU Delft, 2019
%  
%
%  Code used in Manuscript: An automated fracture trace detection technique 
%  using the complex shearlet transform, submitted to Solid Earth
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%  This code extracts fracture traces from images of fractured rock. The
%  extraction is performed using the complex shearlet ridge measure. The
%  ridge measure is computed on an ensemble of shearlet parameters to
%  obtain a ridge ensemble. A highly probable ridge realization is obtained
%  from the ridge ensemble using a threshold ridge strength. The ridges
%  are converted to fractures using image processing steps such as Otsu
%  thresholding, skeletonization and polyline fitting. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%  
%  The code uses functions from the following toolboxes:
%  1. CoSHREM Toolbox by Rafael Reisenhofer 
%       -> CSHRMgetContRidgeSystem ()
%       -> CSHRMsheardec()
%       -> CSHRMgetContRidgeSystem ()

%  2. Geom2D Toolbox by David Legland -> polyline fitting
%  3. MATLAB Mapping Toolbox -> for shapefile I/O operations
%  4. MATLAB Image Processing Toolbox 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
format compact


%% Specify Complex Shearlet Input Parameters    
% Inputs for construction of Shearlet System
% Construct a system of complex-valued shearlets to detect ridges in a 2D grayscale image
% Note: Run this section if you want to create more shearlet systems. 

% ############# LIST OF PARAMETERS FOR SHEARLET CONSTRUCTION #############
% setting the range of shearlet construction parameters and generating combination matrix                                    
waveletEffSuppFraction_Combs = [5 8 12 15];     %   8 - suggested default                               
scalesPerOctave_Combs = [ 1 2 3 ];                %   2 - suggested default 
shearLevel_Combs = [2 3 4];                       %   3 - suggested default   
alpha_Combs = [0 0.5 1];                          % 0.5 - suggested default 
octaves = 3.5;                                    % 3.5 - standard choice 
tic
[ca, cb, cc, cd] = ndgrid(waveletEffSuppFraction_Combs,scalesPerOctave_Combs,shearLevel_Combs,alpha_Combs);
Shearlet_Combs =[ca(:),cb(:),cc(:),cd(:)];
disp(['There are ',num2str(length(Shearlet_Combs)),' combinations of shearlet systems for the given parameter range']);
toc
length_Shearlet_Combs = length(Shearlet_Combs);
clearvars ca cb cc cd 

%% computing all combinations of shearlet systems and saving them in one folder

% specify image dimensions (images are always regular 2D matrices) and the
% resulting shearlet systems are also built for a specific image size
% shearlet systems can be quite large for large image sizes and you may
% run out of memory quite easily. Hence it is advisable to keep the image
% sizes less than 1000 x 1000 pixels. 
rows = 1000; 
cols = 1000;

% specify folder where shearlet systems are saved 
outfolder = 'D:\Github_Test\';

% loop to compute all combinations of shearlet systems
tic
for i=1:3%length(Shearlet_Combs)
    tic    
    waveletEffSupp = ceil(rows/Shearlet_Combs(i,1));
    gaussianEffSupp = ceil(waveletEffSupp/2);
    scalesPerOctave = Shearlet_Combs(i,2);
    shearLevel = Shearlet_Combs(i,3);
    alpha = Shearlet_Combs(i,4);
    scales = scalesPerOctave*octaves;
    %shearletSystem = CSHRMgetContRidgeSystem(rows,cols,ceil(rows/Shearlet_Combs(i,1)),ceil(ceil(rows/Shearlet_Combs(i,1))/2),Shearlet_Combs(i,2),Shearlet_Combs(i,3),Shearlet_Combs(i,4));
    shearletSystem = CSHRMgetContRidgeSystem(rows,cols,waveletEffSupp,gaussianEffSupp,scalesPerOctave,shearLevel,alpha,scales);
    save(strcat([outfolder,'shearletSystem',num2str(i),'.mat']),'shearletSystem','-v7.3');
    disp(['Computing complex shearlet system based on the inputs for combination number: ', num2str(i)]);
    clearvars waveletEffSupp gaussianEffSupp scalesPerOctave shearLevel alpha scales shearletSystem
    toc
end    
toc

%% Setting the range for ridge extraction parameters
% the ridge detection depends on both shearlet system parameter variation
% (which we achieve by using multiple shearlet systems) and ridge parameter
% variations. The total number of ridge realizations is the product of
% 

onlyPositiveOrNegativeRidges = -1;  % +1 -> positive ridges
                                    % -1 -> negative ridges

minContrast_Combs = [5 10 15 20 25];
offset_Combs = [0.001 0.01 0.1 1 2];
[ce, cf] = ndgrid(minContrast_Combs, offset_Combs);
RidgeExtract_Combs =[ce(:),cf(:)];

% number of shearlet systems
%no_Shearlets = length(Shearlet_Combs); 

% if a number of shearlets have already been chosen
no_Shearlets =70;

% total number or ridge realizations 
total_combs = length(RidgeExtract_Combs) * no_Shearlets;
disp(['The total number of ridge realizations is: ', num2str(total_combs)]);
clearvars ce cf offset_Combs minContrast_Combs

%% select fractured rocks image data 
   
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
    
    disp(['You have selected ', num2str(length(InFileList)) , ' images']);

%% This loop computes the ridge realizations and sums them over into one
% normalized ridge image. 

m = 1000;  % rows specified when the shearlet system was built
n = 1000;  % columns specified when the shearlet system was built
    % Note: If the dimensions of the image to be processed are not the same as
    % that of the rows and cols specified when shearlet system was built, the 
    % code doesn't work. The 70 shearlet systems that were used were set for
    % 1000 x 1000 pixel images. If the dimensions don't match, either resize 
    % the raw images (easier) or generate shearlets corresponding to raw image
    % size (takes more time if sizes of all raw images are different)

num_images = length(InFileList); % number of images that are selected

% specify folder where shearlet systems are saved 
outfolder = 'D:\Github_Test\';

% Specify folder where the ridge ensembles are to be saved
output_folder = 'D:\Github_Test\P_Ridges\';

% Loop runs over each image
 for k=1:length(InFileList)
     
  counter = 0;
  tic   
  imageIN = imread(InFileList{k});
  
  % if images have geotiff tags then first 3 channels are extracted 
  % this statement can be commented if all images are grayscale
  %imageIN = rgb2gray(imageIN(:,:,1:3));  
  
  % converting image to grayscale
  imageIN = rgb2gray(imageIN);
  
  % storing the original size of the image     
  [imageIN_Xpixels,imageIN_Ypixels]=size(imageIN);
  
  outfilename = InFileList{k};
  outfilename = InFileListShort{k};
  
  % resizing the image to the rows and columns of shearlet systems
  imageIN = imresize(imageIN,[m n]);
  C_Edges = zeros(m,n,'double');
  
 for i=1:no_Shearlets % number of shearlet systems
     
    tic
     % loading shearlet systems from the shearlet folder
     load(strcat(outfolder,'shearletSystem',num2str(i),'.mat'));
     disp(['Loaded Shearlet System Number: ',num2str(i)]);
    toc
    
    tic
     % creating coefficient matrix for each shearlet (CoSHREM function)
     coeffs = CSHRMsheardec(imageIN,shearletSystem);
     disp(['Created Coeffs Matrix for Shearlet ',num2str(i)]);
    toc
    
  for j=1:15% number of ridge realizations per shearlet
            
    
    tic
    % computing ridge image (CoSHREM function)
    [ridges,~] = CSHRMgetRidges(imageIN,shearletSystem,RidgeExtract_Combs(j,1),RidgeExtract_Combs(j,2),onlyPositiveOrNegativeRidges,coeffs);
    toc
    
    counter = counter + 1;
    
    disp(['Image: ',num2str(k),' ,Shearlet System:', num2str(i),'; Ridge Combination: ',  num2str(counter), ' of ',num2str(total_combs) ]);   

    % In order to check the effects of each ridge realization, you may want
    % to save each ridge to see the effects of each. Can comment these two
    % statements if each intermediate ridge is not required to be saved
    
        %ridges2 = imresize(ridges,[imageIN_Xpixels,imageIN_Ypixels]);
        %imwrite(ridges2, strcat([output_folder,'P_Ridges_',num2str(i),'_','MC',num2str(RidgeExtract_Combs(j,1)),'_Off_',num2str(RidgeExtract_Combs(j,2)),InFileListShort{k}]));


    % summation of computed ridge realization
    C_Edges = C_Edges + ridges;
     
    toc
     clearvars imageOUT ridges
     
  end  
 clearvars shearletSystem coeffs 
 end
 
 C_Edges_norm = mat2gray(C_Edges);
 
 % the summed up ridge is resized to the original image dimension
 C_Edges_norm = imresize(C_Edges_norm,[imageIN_Xpixels,imageIN_Ypixels]);
 
 % writing the ridge ensemble to folder
 imwrite(C_Edges_norm,strcat(output_folder,'P_Ridges_',outfilename));
 disp(['Writing Probabilistic Ridges Map for Image: ',num2str(k), ' of ', num2str(num_images)]);
 clearvars C_Edges_norm C_Edges imageIN outfilename ridges2
 toc
 
 end