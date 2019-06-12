function puls = yapuls(npuls)
% \manchap
%
% Pulsation vector
%
% \mansecSyntax
%
% puls = yapuls(npuls)
%
% \mansecDescription
%
% \libfun{yapuls} returns a pulsation vector \libvar{puls} of size
% \libvar{npuls} which is the concatenation of two subvectors whose
% elements are respectively in  $[0,\pi)$ and $[-\pi,0)$.  Useful function
% when computing wavelets directly in the Fourier domain.
%
% \mansubsecInputData
%
% \begin{description}
% \item[npuls] [REAL SCALAR]: length of the pulsation vector
% \end{description} 
%
% \mansubsecOutputData
%
% \begin{description}
% \item[puls] [REAL VECTOR]: the pulsation vector
% \end{description} 
%
% \mansecExample
%
% \begin{code}
% >> puls5 = yapuls(5)
% >> puls6 = yapuls(6)
% \end{code}
% prints and returns a 5-length, and a 6-length pulsation vector.
%
% \mansecReference
%
% \mansecSeeAlso
%
% vect /cwt.+d$
%
% \mansecLicense
%
% This file is part of YAW Toolbox (Yet Another Wavelet Toolbox)
% You can get it at
% \url{"http://www.fyma.ucl.ac.be/projects/yawtb"}{"yawtb homepage"} 
%
% $Header: /home/cvs/yawtb/tools/misc/yapuls.m,v 1.2 2002-03-26 08:10:25 ljacques Exp $
%
% Copyright (C) 2001-2002, the YAWTB Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

if nargin ~=1,
  error('Argument Mismatch - The function requires one input argument');
end

npuls_2   = floor((npuls-1)/2);
puls      = 2*pi/npuls*[ 0:npuls_2  (npuls_2-npuls+1):-1 ];


% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
