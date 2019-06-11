function L = labelmatrix(CC)
%LABELMATRIX Create label matrix from BWCONNCOMP structure.
%   L = LABELMATRIX(CC) creates a label matrix from the connected
%   components structure CC returned by BWCONNCOMP. 
%
%   The size of L is CC.ImageSize.  The elements of L are integer values
%   greater than or equal to 0.  The pixels labeled 0 are the background.
%   The pixels labeled 1 make up one object, the pixels labeled 2 make up a
%   second object, and so on.
%
%   The class of L depends on CC.NumObjects, and is determined using the
%   following table.
%
%       Class         Range
%       --------      --------------------------------
%       'uint8'       CC.NumObjects <= 255
%       'uint16'      256 <= CC.NumObjects <= 65535
%       'uint32'      65536 <= CC.NumObjects <= 2^32-1
%       'double'      2^32 <= CC.NumObjects 
%
%   LABELMATRIX is more memory-efficient than BWLABEL and BWLABELN because
%   it returns its label matrix in the smallest numeric class necessary for
%   the number of objects.
%
%   Class Support
%   -------------
%   CC is a structure returned by BWCONNCOMP.  L is uint8, uint16, uint32,
%   or double.
%
%   Example
%   -------
%       % Calculate the connected components and display results. Note that
%       % the output from LABELMATRIX uses less memory than the output from
%       % BWLABEL because LABELMATRIX stores its result in the smallest
%       % numeric class necessary for the number of objects.
%
%       BW = imread('text.png');
%       CC = bwconncomp(BW);
%       L = labelmatrix(CC);
%       L2 = bwlabel(BW);
%       whos L L2
%       figure, imshow(label2rgb(L));
%
%   See also BWCONNCOMP, BWLABEL, BWLABELN, LABEL2RGB, REGIONPROPS.

%   Copyright 2008 The MathWorks, Inc.

checkCC(CC, mfilename);

%  'uint8', 'uint16', and 'uint32' are the only integer types labelmatrix
%  can be if we want to maintain memory efficiency when calling label2rgb
%  on the output of labelmatrix.

if CC.NumObjects <= intmax('uint8')
    dataType = 'uint8';
elseif CC.NumObjects <= intmax('uint16')
    dataType = 'uint16';
elseif CC.NumObjects <= intmax('uint32')
    dataType = 'uint32';
else
    dataType = 'double';
end
L = zeros(CC.ImageSize,dataType);

for k = 1 : CC.NumObjects
    L(CC.PixelIdxList{k}) = k;
end