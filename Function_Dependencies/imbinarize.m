function BW = imbinarize(I, varargin)
%IMBINARIZE Binarize image by thresholding.
%   BW = IMBINARIZE(I) binarizes image I with a global threshold computed
%   using Otsu's method, which chooses the threshold to minimize the
%   intraclass variance of the thresholded black and white pixels. BW is
%   the output binary image.
%
%   BW = IMBINARIZE(I, METHOD) binarizes image I with the threshold method
%   specified using METHOD. Available methods are (names can be
%   abbreviated):
%
%   'global'    - Global image threshold using Otsu's method, chosen to
%                 minimize the intraclass variance of the thresholded black
%                 and white pixels. See GRAYTHRESH for details.
%
%   'adaptive'  - Locally adaptive image threshold chosen using local
%                 first-order image statistics around each pixel. See
%                 ADAPTTHRESH for details.
%
%   BW = IMBINARIZE(I, 'adaptive', PARAM1,VAL1,PARAM2,VAL2,...) binarizes
%   image I using name-value pairs to control aspects of adaptive
%   thresholding.
%
%   Parameters include:
%
%   'Sensitivity'           - Specifies the sensitivity factor in the range
%                             [0 1] for adaptive thresholding. A high
%                             sensitivity value leads to thresholding more
%                             pixels as foreground, at the risk of
%                             including some background pixels. Default
%                             value: 0.50
%
%   'ForegroundPolarity'    - Specifies the polarity of the foreground with
%                             respect to the background. Available options
%                             are:
%
%           'bright'        : The foreground is brighter than the
%                             background. (Default)
%           'dark'          : The foreground is darker than the background.
%
%   BW = IMBINARIZE(I, T) binarizes image I using threshold T. T can be a
%   global image threshold, specified as a scalar luminance value or a
%   locally adaptive threshold specified as a matrix of luminance values. T
%   must have values between 0 and 1. If T is a matrix, it must be of the
%   same size as image I. Note that the functions GRAYTHRESH, OTSUTHRESH
%   and ADAPTTHRESH can be used to compute T automatically.
%
%   Class Support
%   -------------
%   The input image I can be a real, non-sparse, 2-D matrix of one of the
%   following classes: uint8, int8, uint16, int16, uint32, int32, single or
%   double. The output binary image BW is a logical matrix of the same size
%   as I.
%
%   Notes
%   -----
%   1. The 'global' method uses a 256-bin image histogram to compute Otsu's
%      threshold. To use a different histogram, see OTSUTHRESH.
%
%   2. The 'adaptive' method binarizes the image using a locally adaptive
%      threshold. A threshold is computed for each pixel using the local
%      mean intensity around the neighborhood of the pixel. This technique
%      is also called Bradley's method. To use a different first order
%      local statistic, see ADAPTTHRESH.
%
%   3. The 'adaptive' method uses a neighborhood size of approximately
%      1/8th of the size of the image (computed as 2*floor(size(I)/16)+1).
%      To use a different neighborhood size, see ADAPTTHRESH.
%
%   4. To produce a binary image from an indexed image, use IND2GRAY to
%      first convert the image to an intensity image. To produce a binary
%      image from an RGB image, use RGB2GRAY to first convert the image to
%      a grayscale intensity image.
%
%   5. IMBINARIZE expects floating point images to be normalized in the
%      range [0,1].
%
%   6. If the image contains Infs or NaNs, the behavior of IMBINARIZE for
%      the 'adaptive' method is undefined. Propagation of Infs or NaNs may
%      not be localized to the neighborhood around Inf/NaN pixels. 
%
%
%   Example 1
%   ---------
%   This example binarizes an image using a global image threshold.
%
%   I = imread('coins.png');
%   BW = imbinarize(I);
%   figure, imshow(BW)
%
%   Example 2
%   ---------
%   This example binarizes an image using a locally adaptive threshold.
%
%   I = imread('rice.png');
%   BW = imbinarize(I, 'adaptive');
%   figure, imshow(BW)
%
%   Example 3
%   ---------
%   This example shows image binarization for images that have darker
%   foreground than background.
%
%   I = imread('printedtext.png');
%   BW = imbinarize(I,'adaptive','ForegroundPolarity','dark','Sensitivity',0.4);
%   figure, imshow(BW)
%
%
%   See also GRAYTHRESH, OTSUTHRESH, ADAPTTHRESH, imageSegmenter.

% Copyright 2015 The MathWorks, Inc.

[I,isNumericThreshold,options] = parseInputs(I,varargin{:});

if isNumericThreshold
    T = options.T;
else
    method = options.Method;
    
    if strcmp(method,'global')
        T = computeGlobalThreshold(I);
    else
        sensitivity = options.Sensitivity;
        fgPolarity  = options.ForegroundPolarity;
        
        T = adaptthresh(I,sensitivity,'ForegroundPolarity',fgPolarity);
    end
    
end

% Binarize image using computed threshold
BW = binarize(I,T);

end

function BW = binarize(I,T)

classrange = getrangefromclass(I);

switch class(I)
    case {'uint8','uint16','uint32'}
        BW = I > T*classrange(2);
        
    case {'int8','int16','int32'}
        BW = I > classrange(1) + (classrange(2)-classrange(1))*T;
        
    case {'single','double'}
        BW = I > T;
end

end

function T = computeGlobalThreshold(I)
% Otsu's threshold is used to compute the global threshold. We convert
% floating point images to uint8 prior to computing the image histogram to
% avoid issues with NaN/Inf values in the input data. im2uint8 nicely
% handles these so that we get a clean histogram for otsuthresh. For other
% types, we compute the histogram in the native type, using 256 bins (this
% is the default in imhist).

if isfloat(I)
    I = im2uint8(I);
    T = otsuthresh( imhist(I) );
else
    T = otsuthresh( imhist(I) );
end
end


%--------------------------------------------------------------------------
% Input Parsing
%--------------------------------------------------------------------------
function [I,isNumericThreshold,options] = parseInputs(I, varargin)

% validate image
validateImage(I);

isNumericThreshold = ~isempty(varargin) && ~ischar(varargin{1});

if isNumericThreshold
    
    options.T = validateT(varargin{1},size(I));
    
    if numel(varargin)>1
        error(message('MATLAB:TooManyInputs'))
    end
    
else
    if isempty(varargin)
        options.Method = 'global';
        return;
    end
    
    options.Method = validatestring(varargin{1},{'global','adaptive'},mfilename,'Method',2);
    
    if strcmp(options.Method,'global')
        
        if numel(varargin)>1
            error(message('MATLAB:TooManyInputs'))
        end
    else
        options.Sensitivity = 0.5;
        options.ForegroundPolarity = 'bright';
        
        numPVArgs = numel(varargin)-1;
        if mod(numPVArgs,2)~=0
            error(message('images:validate:invalidNameValue'));
        end
        
        ParamNames = {'Sensitivity','ForegroundPolarity'};
        ValidateFcn = {@validateSensitivity,@validateForegroundPolarity};
        
        for p = 2 : 2 : numel(varargin)-1
            
            Name = varargin{p};
            Value = varargin{p+1};
            
            idx = strncmpi(Name, ParamNames, numel(Name));
            
            if ~any(idx)
                error(message('images:validate:unknownParamName', Name));
            elseif numel(find(idx))>1
                error(message('images:validate:ambiguousParamName', Name));
            end
            
            validate = ValidateFcn{idx};
            options.(ParamNames{idx}) = validate(Value);
            
        end
    end
    
    
end
end

function validateImage(I)

supportedClasses = {'uint8','uint16','uint32','int8','int16','int32','single','double'};
supportedAttribs = {'real','nonsparse','2d'};
validateattributes(I,supportedClasses,supportedAttribs,mfilename,'I');

end

function T = validateT(T,sizeI)

validateattributes(T,{'numeric'},{'real','nonsparse','2d'},mfilename,'Threshold',2);

if ~( isscalar(T) || isequal(size(T),sizeI) )
    error(message('images:imbinarize:badSizedThreshold'))
end

end

function s = validateSensitivity(s)
validateattributes(s,{'numeric'},{'real','nonsparse','scalar','nonnegative','<=',1},mfilename,'Sensitivity');
end

function fgp = validateForegroundPolarity(fgp)
fgp = validatestring(fgp,{'bright','dark'},mfilename,'ForegroundPolarity');
end
