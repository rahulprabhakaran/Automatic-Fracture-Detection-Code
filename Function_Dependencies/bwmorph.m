function [bwout,lut] = bwmorph(bwin,op,n)
%BWMORPH Morphological operations on binary image.
%   BW2 = BWMORPH(BW1,OPERATION) applies a specific
%   morphological operation to the binary image BW1.
%
%   BW2 = BWMORPH(BW1,OPERATION,N) applies the operation N
%   times.  N can be Inf, in which case the operation is repeated
%   until the image no longer changes.
%
%   OPERATION is a string that can have one of these values:
%      'bothat'       Subtract the input image from its closing
%      'branchpoints' Find branch points of skeleton
%      'bridge'       Bridge previously unconnected pixels
%      'clean'        Remove isolated pixels (1's surrounded by 0's)
%      'close'        Perform binary closure (dilation followed by
%                       erosion)
%      'diag'         Diagonal fill to eliminate 8-connectivity of
%                       background
%      'endpoints'    Find end points of skeleton
%      'fill'         Fill isolated interior pixels (0's surrounded by
%                       1's)
%      'hbreak'       Remove H-connected pixels
%      'majority'     Set a pixel to 1 if five or more pixels in its
%                       3-by-3 neighborhood are 1's
%      'open'         Perform binary opening (erosion followed by
%                       dilation)
%      'remove'       Set a pixel to 0 if its 4-connected neighbors
%                       are all 1's, thus leaving only boundary
%                       pixels
%      'shrink'       With N = Inf, shrink objects to points; shrink
%                       objects with holes to connected rings
%      'skel'         With N = Inf, remove pixels on the boundaries
%                       of objects without allowing objects to break
%                       apart
%      'spur'         Remove end points of lines without removing
%                       small objects completely
%      'thicken'      With N = Inf, thicken objects by adding pixels
%                       to the exterior of objects without connected
%                       previously unconnected objects
%      'thin'         With N = Inf, remove pixels so that an object
%                       without holes shrinks to a minimally
%                       connected stroke, and an object with holes
%                       shrinks to a ring halfway between the hole
%                       and outer boundary
%      'tophat'       Subtract the opening from the input image
%
%   Class Support
%   -------------
%   The input image BW1 can be numeric or logical.
%   It must be 2-D, real and nonsparse.  The output image
%   BW2 is logical.
%
%   Remarks
%   -------
%   To perform erosion or dilation using the structuring element ones(3),
%   use IMERODE or IMDILATE.
%
%   Examples
%   --------
%       BW1 = imread('circles.png');
%       figure, imshow(BW1)
%       BW2 = bwmorph(BW1,'remove');
%       BW3 = bwmorph(BW1,'skel',Inf);
%       figure, imshow(BW2)
%       figure, imshow(BW3)
%
%   See also IMERODE, IMDILATE, BWEULER, BWPERIM.

%   Copyright 1993-2015 The MathWorks, Inc.

% The second output argument, LUT, is intentionally undocumented. In the
% initial release of the Image Processing Toolbox, all the operations
% supported by bwmorph used a single look-up table, which was returned as
% the second output argument.  In subsequent releases, however, bug fixes
% and enhancements resulted in some operations no longer using a single
% look-up table. As a result, the second output argument no longer served
% the purpose envisioned in the original design of the bwmorph syntax. To
% reduce compatibility problems, the second output argument was retained
% in the code, but it has been removed from the documentation.  For
% operations which do not use a single look-up table, the second output
% argument is returned as [].

%
% Input argument parsing
%
if (nargin < 3)
    n = 1;
end

validateattributes(bwin,{'numeric' 'logical'},{'real' 'nonsparse' '2d'}, ...
    mfilename, 'BW', 1);

if ~islogical(bwin)
    bwin = (bwin ~= 0);
end

if ischar(op),
    % BWMORPH(A, 'op', n)
    
    %
    % Find out what operation has been requested
    %
    validOperations = {'bothat',...
        'branchpoints',...
        'bridge',...
        'clean',...
        'close',...
        'diag',...
        'dilate',...
        'endpoints',...
        'erode',...
        'fatten',...
        'fill',...
        'hbreak',...
        'majority',...
        'perim4',...
        'perim8',...
        'open',...
        'remove',...
        'shrink',...
        'skeleton',...
        'spur',...
        'thicken',...
        'thin',...
        'tophat'};
    
    op = validatestring(op,validOperations, 'bwmorph');
    %
    % Call the core function
    %
    
    [bw, lut] = images.internal.algbwmorph(bwin,op,n);    
    
else
    % BWMORPH(A, lut, n)
    
    %
    % Pass on the call to applylut
    %
    lut = op;
    if (isempty(lut))
        error(message('images:bwmorph:emptyLUT'));
    end
    bw   = bwin;
    done = (n <= 0);
    iter = 1;
    while (~done)
        lastbw = bw;
        bw     = applylut(bw, lut);
        done   = ((iter >= n) | isequal(lastbw, bw));
        iter   = iter + 1;
    end
    
   
end

if (nargout == 0)
    imshow(bw);
else
    bwout = bw;
end

if(nargout==2)
    lut = double(lut);
end
if ((nargout == 2) && isempty(lut))
    warning(message('images:bwmorph:lutOutput', op));
end
