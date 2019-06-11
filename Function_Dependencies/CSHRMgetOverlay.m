function imgRgb =  CSHRMgetOverlay(img,edges,forPrint)
% CSHRMgetOverlay Create an overlay image, visualizing the complex shearlet-based ridge or edge
%                 measure on top of the analyzed image.
% 
% Usage (optional parameters are enclosed in angle brackets):
% 
%  imgRgb =  CSHRMgetOverlay(img,edges,<forPrint>)
% 
% Example:
% 
%  img = double(imread('lena.jpg'));
%  shearletSystem = CSHRMgetContEdgeSystem(size(img,1),size(img,2));
%  [edges,orientations] = CSHRMgetEdges(img,shearletSystem);
%  imshow(CSHRMgetOverlay(img,edges));
% 
% 
% See also: CSHRMgetEdges, CSHRMgetRidges
    if nargin < 3
        forPrint = 0;
    end
    if max(edges(:)) > 0
        edges = double(edges)/max(edges(:));
    end
    if forPrint
        imgRgb = 1-(min(0.8,mat2gray(img)*0.4+edges*0.8));
        imgRgb = cat(3,zeros(size(edges)),edges,imgRgb);
        imgRgb = hsl2rgb(imgRgb);
    else
        red = cat(3,edges,zeros(size(edges)),zeros(size(edges)));
        imgRgb = mat2gray(img).*(1-edges);
        imgRgb = cat(3,imgRgb,imgRgb,imgRgb);
        imgRgb = imgRgb + red;
    end    
end

%  Copyright (c) 2016. Rafael Reisenhofer
%
%  Part of CoShREM Toolbox v1.1
%  Built Mon, 11/01/2016
%  This is Copyrighted Material