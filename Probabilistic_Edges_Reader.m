%%
clc
clear all
close all
format compact

%% select probabilistic edges images to convert to binary imaages

    disp('Select multiple image files (they must be in a single folder)');
    [InFileListShort, pathname] = uigetfile('*.jpg;*.tif;*.png','Select input image files','MultiSelect','on');
    
    if not(iscell(InFileListShort))
       InFileListShort = {InFileListShort};
    end
    
    replicatePath = repmat(cellstr(pathname),size(InFileListShort));
    InFileList = strcat(replicatePath,InFileListShort);    
           
    %clear InFileListShort pathname replicatePath;
    funPath = fileparts(which('W125_G63_SPO2_SL3_AL0.5_OCT_3.5_MC00.png'));
    addpath(genpath(funPath));  
    
%%
m=1000;
n=1000;
C_Edges = zeros(m,n,'double');

for i=1:length(InFileList)
    imageIN = imread(InFileList{i});
    imageIN = double(imageIN);
    C_Edges = C_Edges + imageIN;
    i
    clearvars imageIN
end    

%%
C_Edges_norm = mat2gray(C_Edges);
figure(1)
imshow(C_Edges)

figure(2)
imshow(C_Edges_norm)


% using sigmoid nonlinearity
[m,n]=size(C_Edges_norm);
for i=1:m
  for j=1:n  
    if   C_Edges_norm(i,j)~=0
      C_Edges_norm_sigmoid (i,j) = 1 / (1 + exp((-1)*C_Edges_norm(i,j)));
    end  
  end
  i
end
figure(4)
imshow(C_Edges_norm_sigmoid)

%C_Edges_norm_thresh = C_Edges_norm;
C_Edges_norm_thresh = C_Edges_norm_sigmoid;
threshold=0.56;
C_Edges_norm_thresh(C_Edges_norm_thresh<threshold)=0;
C_Edges_norm_thresh(C_Edges_norm_thresh>threshold)=1;
C_Edges_norm_thresh= im2uint8(C_Edges_norm_thresh);
figure(3)
imshow(C_Edges_norm_thresh)
title(['Probability Threshold: ',num2str(threshold)])

%%
% simple linear transformation by mat2gray function
% C_Edges_norm = mat2gray(C_Edges);
% figure(1)
% imshow(C_Edges_norm)
% %clearvars C_Edges_norm
% 
% % using sigmoid nonlinearity
% [m,n]=size(C_Edges_norm);
% for i=1:m
%   for j=1:n  
%     if   C_Edges_norm(i,j)~=0
%       C_Edges_norm_sigmoid (i,j) = 1 / (1 + exp((-1)*C_Edges_norm(i,j)));
%     end  
%   end
%   i
% end
% figure(4)
% imshow(C_Edges_norm_sigmoid)
% 
% % using logit nonlinearity (or reverse sigmoid)
% for i=1:1000
%   for j=1:1000 
%     if   C_Edges_norm(i,j)~=0
%      C_Edges_norm_logit (i,j) = log(C_Edges_norm(i,j)) - log (1 - C_Edges_norm(i,j)) ;
%     end
%   end
%   i
% end
% 
% C_Edges_norm_logit = C_Edges_norm_logit./11.8705;
% figure(5)
% imshow(C_Edges_norm_logit)




%% 
clearvars C_Edges C_Edges_norm C_Edges_norm_sigmoid C_Edges_norm_sigmoid_thresh
output_folder ='D:\PhD\Automatic_Detection\Brejoes_Tif\Brejoes_Binary_Ridges\';
m=1000;
n=1000;
for k=1:length(InFileList)
   tic
   C_Edges = imread(InFileList{k}) ;
   C_Edges_norm = mat2gray(C_Edges);
   for i=1:m
    for j=1:n  
        if   C_Edges_norm(i,j)~=0
          C_Edges_norm_sigmoid (i,j) = 1 / (1 + exp((-1)*C_Edges_norm(i,j)));
        end  
    end
    i
   end   
    threshold=0.595;
    C_Edges_norm_sigmoid_thresh = C_Edges_norm_sigmoid;
    C_Edges_norm_sigmoid_thresh(C_Edges_norm_sigmoid_thresh<threshold)=0;
    C_Edges_norm_sigmoid_thresh(C_Edges_norm_sigmoid_thresh>threshold)=1;
    % C_Edges_norm_sigmoid_thresh= im2uint8(C_Edges_norm_sigmoid_thresh);
    % colormap hot
    figure(9)
    imshow(C_Edges_norm_sigmoid_thresh)
    % set(gca,'DefaultTextFontSize',20)
    % title(['Intensity Threshold: ',num2str(threshold)])
    output_filename = InFileListShort{k};
    output_filename = strcat('B_Ridges_',output_filename(9:21));
    imwrite(C_Edges_norm_sigmoid_thresh,strcat(output_folder,output_filename));
    k
    toc
end    
%InFileList{k}


%figure(4)
%imshow(C_Edges_norm_sigmoid)
%colormap hot

%%


%%
%imwrite(C_Edges_norm_sigmoid_thresh,'D:\PhD\Automatic_Detection\Brejoes_Tif\Brejoes_Binary_Ridges\B_Ridges_0006_0030.tif'); 

