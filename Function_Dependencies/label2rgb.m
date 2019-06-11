function RGB = label2rgb(varargin) 
%LABEL2RGB Convert label matrix to RGB image.
%   RGB = LABEL2RGB(L) converts a label matrix L, such as returned by
%   LABELMATRIX, BWLABEL, BWLABELN, or WATERSHED, into a color RGB image
%   for the purpose of visualizing the labeled regions.
%
%   RGB = LABEL2RGB(L, MAP) defines the colormap to be used in the RGB
%   image.  MAP can either be an n x 3 colormap matrix, a string containing
%   the name of a colormap function (such as 'jet' or 'gray'), or a
%   function handle of a colormap function (such as @jet or @gray).
%   LABEL2RGB evaluates MAP so that there is a different color for each
%   region in L. If MAP is not specified, 'jet' is used as the default.
%
%   RGB = LABEL2RGB(L, MAP, ZEROCOLOR) defines the RGB color of the
%   elements labeled 0 in the input label matrix L.  ZEROCOLOR can either
%   be an RGB triple, or one of the following: 'y' (yellow), 'm',
%   (magenta), 'c' (cyan), 'r'(red), 'g' (green), 'b' (blue), 'w' (white),
%   or 'k' (black). If ZEROCOLOR is not specified, [1 1 1] is used as the
%   default.
%   
%   RGB = LABEL2RGB(L, MAP, ZEROCOLOR, ORDER), controls how colormap colors
%   are assigned to regions in the label matrix.  If ORDER is 'noshuffle'
%   (the default), then colormap colors are assigned to the label matrix
%   regions in numerical order.  If ORDER is 'shuffle', then colormap
%   colors are pseudorandomly shuffled.
%   
%   Class Support
%   -------------
%   The input label matrix L can have any numeric class. It must contain
%   finite nonnegative integers.  RGB is uint8.
%
%   Example 1
%   ---------
%   %Use label2rgb to customize display of label matrix.
%
%       I = imread('rice.png');
%       figure, imshow(I)
%       BW = imbinarize(I);
%       CC = bwconncomp(BW);
%       L = labelmatrix(CC);
%       RGB = label2rgb(L);
%       RGB2 = label2rgb(L, 'spring', 'c', 'shuffle');
%       figure, imshow(RGB), figure, imshow(RGB2)
%
%   See also BWCONNCOMP,BWLABEL,COLORMAP,ISMEMBER,LABELMATRIX,WATERSHED.

%   Copyright 1993-2015 The MathWorks, Inc.

[label,map,zerocolor,order,fcnflag] = parse_inputs(varargin{:});

% Determine the number of regions in the label matrix.
numregion = double(max(label(:)));

% If MAP is a function, evaluate it.  Make sure that the evaluated function
% returns a valid colormap.
if  fcnflag == 1
    if numregion == 0
      cmap = [];
    else
      cmap = feval(map, numregion);
      if ~isreal(cmap) || any(cmap(:) > 1) || any(cmap(:) < 0) || ...
            ~isequal(size(cmap,2),3) || size(cmap,1) < 1
        error(message('images:label2rgb:functionReturnsInvalidColormap'));
      end
    end
else
    cmap = map;
end

% If ORDER is set to 'shuffle', create a private stream with a fixed seed,
% which creates the same "random" permutation every time it is called.
if isequal(order,'shuffle')
    stream = RandStream('swb2712','seed',0);
    index = randperm(stream,numregion);
    cmap = cmap(index,:,:);
end

% Issue a warning if the zerocolor (boundary color) matches the color of one
% of the regions. 
for i=1:numregion
  if isequal(zerocolor,cmap(i,:))
    warning(message('images:label2rgb:zerocolorSameAsRegionColor', i));
  end
end
cmap = [zerocolor;cmap];

if isa(label,'uint8') || isa(label,'uint16') || isa(label,'uint32')
    RGB = ind2rgb8(label, cmap);
else
    % Using label + 1 for two reasons: 1) IND2RGB and IND2RGB8 do not like
    % double arrays containing zero values, and 2)for non-double, IND2RGB would
    % cast to a double and do this.
    RGB = ind2rgb8(double(label)+1,cmap);
end
 
%  Function: parse_inputs
%  ----------------------
function [L, Map, Zerocolor, Order, Fcnflag] = parse_inputs(varargin) 
% L         label matrix: matrix containing non-negative values.  
% Map       colormap: name of standard colormap, user-defined map, function
%           handle.
% Zerocolor RGB triple or Colorspec
% Order     keyword if specified: 'shuffle' or 'noshuffle'
% Fcnflag   flag to indicating that Map is a function


narginchk(1,4);

% set defaults
L = varargin{1};
Map = 'jet';    
Zerocolor = [1 1 1]; 
Order = 'noshuffle';
Fcnflag = 0;

% parse inputs
if nargin > 1
    Map = varargin{2};
end
if nargin > 2
    Zerocolor = varargin{3};
end
if nargin > 3
    Order = varargin{4};
end

% error checking for L
validateattributes(L,{'numeric','logical'}, ...
              {'real' '2d' 'nonsparse' 'finite' 'nonnegative' 'integer'}, ...
              mfilename,'L',1);

% error checking for Map
[fcn, fcnchk_msg] = fcnchk(Map);
if isempty(fcnchk_msg)
    Map = fcn;
    Fcnflag = 1;
else
    if isnumeric(Map)
        if ~isreal(Map) || any(Map(:) > 1) || any(Map(:) < 0) || ...
                    ~isequal(size(Map,2), 3) || size(Map,1) < 1
          error(message('images:label2rgb:invalidColormap'));
        end
    else
        error(fcnchk_msg);
    end
end    
    
% error checking for Zerocolor
if ~ischar(Zerocolor)
    % check if Zerocolor is a RGB triple
    if ~isreal(Zerocolor) || ~isequal(size(Zerocolor),[1 3]) || ...
                any(Zerocolor> 1) || any(Zerocolor < 0)
      error(message('images:label2rgb:invalidZerocolor'));
    end
else    
    [cspec, msg] = cspecchk(Zerocolor);
    if ~isempty(msg)
	%message is translated at source.
        error(message('images:label2rgb:notInColorspec', msg))
    else
        Zerocolor = cspec;
    end
end

% error checking for Order
valid_order = {'shuffle', 'noshuffle'};
idx = strncmpi(Order, valid_order,length(Order));
if ~any(idx)
    error(message('images:label2rgb:invalidEntryForOrder'))
elseif nnz(idx) > 1
    error(message('images:label2rgb:ambiguousEntryForOrder', Order))
else
    Order = valid_order{idx};
end
