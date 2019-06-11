function CC = bwconncomp(varargin)
%BWCONNCOMP Find connected components in binary image.
%   CC = BWCONNCOMP(BW) returns the connected components CC found in BW. 
%   BW is a binary image that can have any dimension. CC is a structure 
%   with four fields:
%
%      Connectivity   Connectivity of the connected components (objects).
%
%      ImageSize      Size of BW.
%
%      NumObjects     Number of connected components (objects) in BW.
%
%      PixelIdxList   1-by-NumObjects cell array where the kth element
%                     in the cell array is a vector containing the linear
%                     indices of the pixels in the kth object.
%    
%   BWCONNCOMP uses a default connectivity of 8 for two dimensions, 26 for
%   three dimensions, and CONNDEF(NDIMS(BW),'maximal') for higher
%   dimensions.
%
%   CC = BWCONNCOMP(BW,CONN) specifies the desired connectivity for
%   the connected components.  CONN may have the following scalar values:
%
%      4             two-dimensional four-connected neighborhood
%      8             two-dimensional eight-connected neighborhood
%      6             three-dimensional six-connected neighborhood
%      18            three-dimensional 18-connected neighborhood
%      26            three-dimensional 26-connected neighborhood
%     
%   Connectivity may be defined in a more general way for any dimension 
%   using a 3-by-3-by- ... -by-3 matrix of 0s and 1s.  The
%   1-valued elements define neighborhood locations relative to the center
%   element of CONN.  CONN must be symmetric about its center element.
%
%   Note: On the use of BWLABEL, BWLABELN, BWCONNCOMP, and REGIONPROPS
%   ------------------------------------------------------------------
%   The functions BWLABEL, BWLABELN, and BWCONNCOMP all compute connected
%   components for binary images.  BWCONNCOMP is the most recent addition
%   to the Image Processing Toolbox and is intended to replace the use
%   of BWLABEL and BWLABELN.  It uses significantly less memory and is
%   sometimes faster than the older functions.
%  
%                Input  Output            Memory   Connectivity
%                Dim    Form              Use
%                ----------------------------------------------
%   BWLABEL      2-D    Double-precision  High     4 or 8
%                       label matrix           
% 
%   BWLABELN     N-D    Double-precision  High     Any
%                       label matrix       
% 
%   BWCONNCOMP   N-D    CC struct         Low      Any
%
%   To extract features from a binary image using REGIONPROPS using the
%   default connectivity, just pass BW directly into REGIONPROPS, i.e.,
%   REGIONPROPS(BW). 
%
%   To compute a label matrix having more memory-efficient data type (e.g.
%   uint8 versus double), use the LABELMATRIX function on the output of
%   BWCONNCOMP. See the documentation for each function for more
%   information.
%
%   Class Support
%   -------------
%   BW can be logical or numeric, and it must be real, N-D, and nonsparse.
%   CC is a structure.
%
%   Example 1
%   ---------
%   % Calculate the centroids of the 3D objects.
%
%       BW = cat(3,[1 1 0; 0 0 0; 1 0 0],...
%                  [0 1 0; 0 0 0; 0 1 0],...
%                  [0 1 1; 0 0 0; 0 0 1])
%       CC = bwconncomp(BW);
%       S = regionprops(CC,'Centroid');
%
%   Example 2
%   ---------
%   % Erase the biggest object from the image.
%
%       BW = imread('text.png');
%       imshow(BW);
%
%       CC = bwconncomp(BW);
%       numPixels = cellfun(@numel,CC.PixelIdxList);
%       [biggest,idx] = max(numPixels);
%       BW(CC.PixelIdxList{idx}) = 0;
%
%       figure, imshow(BW);  
%
%   See also BWLABEL, BWLABELN, LABELMATRIX, REGIONPROPS.

%   Copyright 2008-2015 The MathWorks, Inc.

[BW, conn] = parseInputs(varargin{:});

CC = struct(...
    'Connectivity', conn, ...
    'ImageSize', size(BW), ...
    'NumObjects', [], ...
    'PixelIdxList', []);

if ismatrix(BW) && (isequal(conn,4) || isequal(conn,8))
    [CC.PixelIdxList,CC.NumObjects] = bwconncomp_2d(BW, conn);
else
    [CC.PixelIdxList,CC.NumObjects] = bwconncomp_nd(BW, conn);
end

%--------------------------------------------------------------------------
function [BW,conn] = parseInputs(varargin)

narginchk(1,2);

BW = varargin{1};
validateattributes(BW, {'logical' 'numeric'}, {'real', 'nonsparse'}, ...
              mfilename, 'BW', 1);
if ~islogical(BW)
    BW = BW ~= 0;
end

if nargin < 2
    %BWCONNCOMP(BW)
    % special case 8 and 26 because these numbers would be more
    % understandable when returned to the user in the connectivity field as
    % opposed to returning ones(3) or ones(3,3,3) for the connectivity.
    % Fits in better with bwlabel and bwlabeln.
    if ismatrix(BW)
        conn = 8;
    elseif ndims(BW) == 3
        conn = 26;
    else
        conn = conndef(ndims(BW), 'maximal');
    end
else
    %BWCONNCOMP(BW,CONN)
    conn = varargin{2};
    iptcheckconn(conn,mfilename,'CONN',2);
    
    % special case so that we go through the 2D code path for 4 or 8
    % connectivity
    if isequal(conn, [0 1 0;1 1 1;0 1 0])
        conn = 4;
    end
    if isequal(conn, ones(3))
        conn = 8;
    end
end
