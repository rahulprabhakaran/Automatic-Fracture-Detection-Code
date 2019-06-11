function ags = Angles2D(lns)
% Angles2D
% returns angles of 2D lines (fracture traces)
% the angles are always returned in counterclockwise direction if y2 > y1
% 
% example Angles2D ([0 0  1  1])  returns   45 deg
%         Angles2D ([0 0 -1  1])  returns  135 deg
%         Angles2D ([0 0  1 -1])  returns  -45 deg
%         Angles2D ([0 0 -1 -1])  returns -135 deg
% 
%         Angles2D ([ 1  1 0 0])  returns -135 deg
%         Angles2D ([-1  1 0 0])  returns  -45 deg
%         Angles2D ([ 1 -1 0 0])  returns  135 deg
%         Angles2D ([-1 -1 0 0])  returns   45 deg
% Usage :
% ags = Angles2D(lns)
%
% input : lns       (n, 4)
% output: ags       (n) radian [-pi..pi]
%
% Alghalandis Discrete Fracture Network Engineering (ADFNE)
% Author: Younes Fadakar Alghalandis
% Copyright (c) 2016 Alghalandis Computing @ http://alghalandis.net
% All rights reserved.

ags = atan2(lns(:, 4)-lns(:, 2), lns(:, 3)-lns(:, 1));

% converting radians to degrees
ags = ags.*57.2958;

