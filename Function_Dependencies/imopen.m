function B = imopen(A, se_) %#codegen
%IMOPEN Morphologically open image.
%   IM2 = IMOPEN(IM,SE) performs morphological opening on the grayscale
%   or binary image IM with the structuring element SE.  SE must be a
%   single structuring element object, as opposed to an array of
%   objects.
%
%   IM2 = IMOPEN(IM,NHOOD) performs opening with the structuring element
%   STREL(NHOOD), where NHOOD is an array of 0s and 1s that specifies the
%   structuring element neighborhood.
%
%   The morphological open operation is an erosion followed by a dilation,
%   using the same structuring element for both operations.
%
%   Class Support
%   -------------
%   IM can be any numeric or logical class and any dimension, and must be
%   nonsparse.  If IM is logical, then SE must be flat.  IM2 has the same
%   class as IM.
%
%   Example
%   -------
%   Remove snowflakes having a radius less than 5 pixels by opening it with
%   a disk-shaped structuring element having a 5 pixel radius.
%
%       original = imread('snowflakes.png');
%       se = strel('disk',5);
%       afterOpening = imopen(original,se);
%       figure, imshow(original), figure, imshow(afterOpening,[])
%
%   See also IMCLOSE, IMDILATE, IMERODE, STREL.

%   Copyright 1993-2015 The MathWorks, Inc.

validateattributes(A, ...
    {'numeric' 'logical'}, ...
    {'real' 'nonsparse'},...
    mfilename, 'I or BW', 1);

se = images.internal.strelcheck(se_, mfilename, 'SE', 2);

coder.internal.errorIf((length(se(:)) > 1),...
    'images:imopen:nonscalarStrel');

strel_is_flat    = coder.const(isflat(se));
input_is_logical = coder.const(islogical(A));
strel_is_2d      = coder.const(ismatrix(getnhood(se)));

coder.internal.errorIf((input_is_logical && ~strel_is_flat), ...
    'images:imopen:binaryImageWithNonflatStrel');

input_is_2d = coder.const(numel(size(A))==2);

pre_pack_ = coder.const(input_is_logical && input_is_2d && strel_is_2d);
coder.extrinsic('images.internal.coder.isCodegenForHost');
coder.extrinsic('images.internal.coder.useSharedLibrary');
if(coder.target('MATLAB'))
    pre_pack = pre_pack_;
else    
    % packed inputs are only supported in shared library mode (host)
    pre_pack = pre_pack_ ...
        && ~coder.const(images.internal.coder.isCodegenForHost())...
        && coder.const(images.internal.coder.useSharedLibrary());
end

pre_pack = coder.const(pre_pack);

M = size(A,1);
if pre_pack    
    inputImage = bwpack(A);
    packopt = 'ispacked';
    outputImage = imdilate(imerode(inputImage,se,packopt,M),se,packopt,M);
    B = bwunpack(outputImage,M);
else
    inputImage = A;
    packopt = 'notpacked';
    B = imdilate(imerode(inputImage,se,packopt,M),se,packopt,M);
end