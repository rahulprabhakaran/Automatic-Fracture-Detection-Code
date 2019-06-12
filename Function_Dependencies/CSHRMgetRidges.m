function [ridges,tangentOrientations,widths,heights] = CSHRMgetRidges(img,shearletSystem,minContrast,offset,onlyPositiveOrNegativeRidges,coeffs)
% CSHRMgetRidges Compute the complex shearlet-based ridge measure of a 2D grayscale image.
% 
% Usage (optional parameters are enclosed in angle brackets):
% 
%  [ridges,tangentOrientations,widths,heights] = CSHRMgetRidges(img,sys,<minContrast>,<offset>,<scalesUsedForPivotSearch>)
% 
% Input:
%
%                           img: 2D grayscale image.
%                           sys: A complex shearlet system. Can be constructed using CSHRMgetContRidgeSystem.
%                   minContrast: (optional) Specifies the minimal contrast between a ridge and its background for the ridge to be detected. Default:                             
%                        offset: (optional) Determines the difference in octaves between the scales of odd-symmetric shearlets and 
%                                even-symmetric shearlets used to compute the complex shearlet-based ridge measure. Default: 1.
%  onlyPositiveOrNegativeRidges: (optional) If set to 1, only positive ridges are detected. If set to -1, only negative ridges are detected.
%  
% 
% Output:
% 
%               ridges: The complex shearlet-based ridge measure computed each pixel of the analyzed image. The values are ranging
%                       from 0 to 1. 
%  tangentOrientations: Approximations of the local tangent orientations with values ranging from 0? to 180?.
%               widths: Width-estimates of the detected ridges in pixels. 
%              heights: Heigh-estimates of the detected ridges.
% Example:
% 
%  img = double(imread('lena.jpg'));
%  shearletSystem = CSHRMgetContRidgeSystem(size(img,1),size(img,2));
%  [ridges,tangentOrientations,widths,heights] = CSHRMgetRidges(img,shearletSystem);
%  imshow(CSHRMgetOverlay(img,ridges));
% 
% 
% See also: CSHRMgetEdges, CSHRMgetContRidgeSystem    
    if (nargin < 3), minContrast = (max(img(:)) - min(img(:)))/30; end
    if (nargin < 4), offset = 1; end
    if (nargin < 5), onlyPositiveOrNegativeRidges = 0; end
    
    %coeffs = CSHRMsheardec(img,shearletSystem);
    offset = floor(offset*shearletSystem.scalesPerOctave);
    nOrientations = 2*(2^shearletSystem.shearLevel+2);
    coeffs = reshape(coeffs,size(coeffs,1),size(coeffs,2),nOrientations,size(coeffs,3)/(nOrientations));

    ci = abs(imag(coeffs(:,:,:,1:(end-offset))));
    cr = real(coeffs(:,:,:,(1+offset):end));
    
    [~,maxPivot] = max(reshape(abs(cr),size(cr,1),size(cr,2),size(cr,3)*size(cr,4)),[],3);
    pivotOris = mod(maxPivot-1,nOrientations)+1;
    pivotScales = fix((maxPivot-1)/nOrientations)+1;    

    positiveWidths = uint16(shearletSystem.positiveWidths((offset*2+1):end));
    expectedCoeffs = shearletSystem.expectedCoeffs(1:max(positiveWidths),(offset+1):end);
    
    [ridges,tangentOrientations,widths,heights] = CSHRMgetRidgesAndTangentOrientationsMex(cr,ci,uint16(pivotOris),uint16(pivotScales),positiveWidths,expectedCoeffs,shearletSystem.shearLevel,minContrast);
    tangentOrientations(tangentOrientations >= 0) = CSHRMmapOrientationsToAngles(tangentOrientations(tangentOrientations >= 0),shearletSystem.shearLevel);

    if onlyPositiveOrNegativeRidges == 1
        ridges = ridges.*(heights > 0);
    elseif onlyPositiveOrNegativeRidges == -1
        ridges = ridges.*(heights < 0);
    end
    
    helpIndexes = ridges<eps;
    tangentOrientations(helpIndexes) = -1;
    widths(helpIndexes) = 0;
    heights(helpIndexes) = 0;
    
end


%  Copyright (c) 2016. Rafael Reisenhofer
%
%  Part of CoShREM Toolbox v1.1
%  Built Mon, 11/01/2016
%  This is Copyrighted Material