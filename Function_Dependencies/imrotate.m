function varargout = imrotate(varargin)
%IMROTATE Rotate image.
%   B = IMROTATE(A,ANGLE) rotates image A by ANGLE degrees in a
%   counterclockwise direction around its center point. To rotate the image
%   clockwise, specify a negative value for ANGLE. IMROTATE makes the output
%   image B large enough to contain the entire rotated image. IMROTATE uses
%   nearest neighbor interpolation, setting the values of pixels in B that
%   are outside the rotated image to 0 (zero).
%
%   B = IMROTATE(A,ANGLE,METHOD) rotates image A, using the interpolation
%   method specified by METHOD. METHOD is a string that can have one of the
%   following values. The default value is enclosed in braces ({}).
%
%        {'nearest'}  Nearest neighbor interpolation
%
%        'bilinear'   Bilinear interpolation
%
%        'bicubic'    Bicubic interpolation. Note: This interpolation
%                     method can produce pixel values outside the original
%                     range.
%
%   B = IMROTATE(A,ANGLE,METHOD,BBOX) rotates image A, where BBOX specifies
%   the size of the output image B. BBOX is a text string that can have
%   either of the following values. The default value is enclosed in braces
%   ({}).
%
%        {'loose'}    Make output image B large enough to contain the
%                     entire rotated image. B is generally larger than A.
%
%        'crop'       Make output image B the same size as the input image
%                     A, cropping the rotated image to fit.
%
%   Class Support
%   -------------
%   The input image can be numeric or logical.  The output image is of the
%   same class as the input image.
%
%   Note
%   ----
%   The function IMROTATE changed in version 9.3 (R2015b).  Previous 
%   versions of the Image Processing Toolbox use different spatial 
%   conventions.  If you need the same results produced by the previous 
%   implementation, use the function IMROTATE_OLD.
%
%   This function may take advantage of hardware optimization for datatypes
%   uint8, uint16, single, and double to run faster.
%
%   Example
%   -------
%        % This example brings image I into horizontal alignment by
%        % rotating the image by -1 degree.
%
%        I = fitsread('solarspectra.fts');
%        I = mat2gray(I);
%        J = imrotate(I,-1,'bilinear','crop');
%        figure, imshow(I), figure, imshow(J)
%
%   See also IMROTATE_OLD, IMCROP, IMRESIZE, IMTRANSFORM, TFORMARRAY.

%   Copyright 1992-2016 The MathWorks, Inc.

% Grandfathered:
%   Without output arguments, IMROTATE(...) displays the rotated
%   image in the current axis.

[A,ang,method,bbox] = parse_inputs(varargin{:});

if (isempty(A))
    
    B = A; % No rotation needed
    
else
    so = size(A);
    twod_size = so(1:2);
    
    if rem(ang,90) == 0
        % Catch and speed up 90 degree rotations
        
        % determine if angle is +- 90 degrees or 0,180 degrees.
        multiple_of_ninety = mod(floor(ang/90), 4);
        
        % initialize array of subscripts
        v = repmat({':'},[1 ndims(A)]);
        
        switch multiple_of_ninety
            
            case 0
                % 0 rotation;
                B = A;
                
            case {1,3}
                % +- 90 deg rotation
                
                thirdD = prod(so(3:end));
                A = reshape(A,[twod_size thirdD]);
                
                not_square = twod_size(1) ~= twod_size(2);
                if strcmpi(bbox, 'crop') && not_square
                    % center rotated image and preserve size
                    
                    imbegin = (max(twod_size) == so)*abs(diff(floor(twod_size/2)));
                    vec = 1:min(twod_size);
                    v(1) = {imbegin(1)+vec};
                    v(2) = {imbegin(2)+vec};
                    
                    new_size = [twod_size thirdD];
                    
                else
                    % don't preserve original size
                    new_size = [fliplr(twod_size) thirdD];
                end
                
                % pre-allocate array
                if islogical(A)
                    B = false(new_size);
                else
                    B = zeros(new_size, 'like', A);
                end
                
                B(v{1},v{2},:) = rot90(A(v{1},v{2},:), multiple_of_ninety);
                
                B = reshape(B,[new_size(1) new_size(2) so(3:end)]);
                
            case 2
                % 180 rotation
                
                v(1) = {twod_size(1):-1:1};
                v(2) = {twod_size(2):-1:1};
                B = A(v{:});
        end
        
    else % Perform general rotation
        
        tform = affine2d([cosd(ang) -sind(ang) 0; sind(ang) cosd(ang) 0; 0 0 1]);
        
        if useIPP(A)
            % The library routine has different edge behavior than our code.
            % This difference can be worked around with zero padding.
            % Zero padding is handled by the IPP symbols
            outputSize = getOutputBound(tform,twod_size,bbox);
            
            if isreal(A)
                B = imrotatemex(A,ang,outputSize,method);
            else
                % If input is complex valued, call imrotatemex on real and
                % imaginary parts separately then combine results.
                B = complex(imrotatemex(real(A),ang,outputSize,method),...
                    imrotatemex(imag(A),ang,outputSize,method));
            end
            
        else % rotate using IMWARP
            
            RA = imref2d(size(A));
            Rout = images.spatialref.internal.applyGeometricTransformToSpatialRef(RA,tform);
            
            if strcmp(bbox,'crop')
                % Trim Rout, preserve center and resolution.
                Rout.ImageSize = RA.ImageSize;
                xTrans = mean(Rout.XWorldLimits) - mean(RA.XWorldLimits);
                yTrans = mean(Rout.YWorldLimits) - mean(RA.YWorldLimits);
                Rout.XWorldLimits = RA.XWorldLimits+xTrans;
                Rout.YWorldLimits = RA.YWorldLimits+yTrans;
            end
            
            B = imwarp(A,tform,method,'OutputView',Rout, 'SmoothEdges',true);
        end
    end
    
end


% Output
switch nargout
    case 0
        % Need to set varargout{1} so ans gets populated even if user doesn't ask for output
        varargout{1} = B;
    case 1
        varargout{1} = B;
    case 3
        error(message('images:removed:syntax','[R,G,B] = IMROTATE(RGB)','RGB2 = IMROTATE(RGB1)'))
    otherwise
        error(message('images:imrotate:tooManyOutputs'))
end
end


function TF = useIPP(A)

% We enable acceleration for uint8, uint16, and single or double precision
% floating point inputs.
supportedType = isa(A,'uint8') || isa(A,'uint16') || isa(A,'float');
% IPP's rotate behavior for vectors is different than the non IPP version
% Defaulting to the non-IPP code path for vector inputs
TF =  images.internal.useIPPLibrary() && supportedType && ~isProblemSizeTooBig(A) && ~isvector(A); 
end

function TF = isProblemSizeTooBig(A)

% IPP cannot handle double-precision inputs that are too big. Switch to
% using tform when the image is double-precision and is too big.

imageIsDoublePrecision = isa(A,'double');

padSize = 2;
numel2DInputImage = (size(A,1) + 2*padSize) * (size(A,2) + 2*padSize);

% The size threshold is double(intmax('int32'))/8. The double-precision
% IPP routine can only handle images that have fewer than this many pixels.
% This is hypothesized to be because they use an int to hold a pointer
% offset for the input image. This overflows when the offset becomes large
% enough that ptrOffset*sizeof(double) exceeds intmax.
sizeThreshold = 2.6844e+08;
TF = imageIsDoublePrecision && (numel2DInputImage>=sizeThreshold);
end

function [A,ang,method,bbox] = parse_inputs(varargin)

narginchk(2,4);

% validate image
A = varargin{1};
validateattributes(A,{'numeric','logical'},{},mfilename,'input image',1);

% validate angle
ang = double(varargin{2});
validateattributes(ang,{'numeric'},{'real','scalar'},mfilename,'ANGLE',2);

method = 'nearest';
bbox   = 'loose';
strings  = {'nearest','bilinear','bicubic','crop','loose'};
isBBox   = [ false   ,false     ,false    ,true  ,true   ];
if nargin==3
    
    arg = varargin{3};
    if ~ischar(arg)
        error(message('images:imrotate:expectedString'));
    end
    idx = stringmatch(lower(arg),strings);
    checkStringValidity(idx,arg);
    arg = strings{idx};
    
    if isBBox(idx)
        bbox = arg;
    else
        method = arg;
    end
    
elseif nargin==4
    
    arg1 = varargin{3};
    if ~ischar(arg1)
        error(message('images:imrotate:expectedString'));
    end
    idx1 = stringmatch(lower(arg1),strings);
    checkStringValidity(idx1,arg1);
    arg1 = strings{idx1};
    
    arg2 = varargin{4};
    if ~ischar(arg2)
        error(message('images:imrotate:expectedString'));
    end
    idx2 = stringmatch(lower(arg2),strings);
    checkStringValidity(idx2,arg2);
    arg2 = strings{idx2};
    
    if isBBox(idx1)
        bbox = arg1;
    else
        method = arg1;
    end
    
    if isBBox(idx2)
        bbox = arg2;
    else
        method = arg2;
    end
end
end

function idx = stringmatch(str,cellOfStrings)
idx = find(strncmpi(str, cellOfStrings, numel(str)));
end

function checkStringValidity(idx,arg)
if isempty(idx)
    error(message('images:imrotate:unrecognizedInputString', arg));
elseif numel(idx)>1
    error(message('images:imrotate:ambiguousInputString', arg));
end
end

function outputSize = getOutputBound(tform,twod_size,bbox)

XWorldLimits = [0.5 twod_size(2)+0.5];
YWorldLimits = [0.5 twod_size(1)+0.5];

if strcmpi(bbox,'loose')
    
    [XWorldLimitsOut,YWorldLimitsOut] = outputLimits(tform,XWorldLimits,YWorldLimits);
    
    % Use ceil to provide grid that will accomodate world limits at roughly
    % the target resolution. Assumes that the resolution is 1.
    numCols = ceil(diff(XWorldLimitsOut));
    numRows = ceil(diff(YWorldLimitsOut));
    
    outputSize = [numRows numCols];
else
    outputSize = twod_size;
end
end