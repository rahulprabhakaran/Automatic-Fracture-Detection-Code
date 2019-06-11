function I = mat2gray(A,limits)
%MAT2GRAY Convert matrix to intensity image.
%   I = MAT2GRAY(A,[AMIN AMAX]) converts the matrix A to the intensity image I.
%   The returned matrix I contains values in the range 0.0 (black) to 1.0 (full
%   intensity or white).  AMIN and AMAX are the values in A that correspond to
%   0.0 and 1.0 in I.  Values less than AMIN become 0.0, and values greater than
%   AMAX become 1.0.
%
%   I = MAT2GRAY(A) sets the values of AMIN and AMAX to the minimum and maximum
%   values in A.
%
%   Class Support
%   -------------  
%   The input array A can be logical or numeric. The output image I is double.
%
%   Example
%   -------
%       I = imread('rice.png');
%       J = filter2(fspecial('sobel'), I);
%       K = mat2gray(J);
%       figure, imshow(I), figure, imshow(K)
%
%   See also GRAY2IND, IND2GRAY, RGB2GRAY.

%   Copyright 1992-2013 The MathWorks, Inc.


validateattributes(A,{'logical','uint8', 'uint16', 'uint32',...
                    'int8', 'int16', 'int32','single', 'double'},...
                    {},mfilename,'A',1);

if nargin == 1
  limits = double([min(A(:)) max(A(:))]);
else
  validateattributes(limits,{'double'},{'numel',2},mfilename,'LIMITS',2);
end

if limits(2)==limits(1)   % Constant Image
   I = double(A);
else
  delta = 1 / (limits(2) -  limits(1));
  I = imlincomb(delta, A, -limits(1)*delta, 'double');
end

% Make sure all values in I are between 0 and 1.
I = max(0,min(I,1));   
