function result = DouglasPeucker(Points,epsilon)
% The Ramer–Douglas–Peucker algorithm (RDP) is an algorithm for reducing 
% the number of points in a curve that is approximated by a series of 
% points. The initial form of the algorithm was independently suggested 
% in 1972 by Urs Ramer and 1973 by David Douglas and Thomas Peucker and 
% several others in the following decade. This algorithm is also known 
% under the names Douglas–Peucker algorithm, iterative end-point fit 
% algorithm and split-and-merge algorithm. [Source Wikipedia]
%
% Input:
%           Points: List of Points 2xN
%           epsilon: distance dimension, specifies the similarity between
%           the original curve and the approximated (smaller the epsilon,
%           the curves more similar)
% Output:
%           result: List of Points for the approximated curve 2xM (M<=N)    
%           
%
% -------------------------------------------------------
% Code: Reza Ahmadzadeh (2017) 
% -------------------------------------------------------
dmax = 0;
edx = length(Points);
for ii = 2:edx-1
    d = penDistance(Points(:,ii),Points(:,1),Points(:,edx));
    if d > dmax
        idx = ii;
        dmax = d;
    end
end

if dmax > epsilon
    % recursive call
    recResult1 = DouglasPeucker(Points(:,1:idx),epsilon);
    recResult2 = DouglasPeucker(Points(:,idx:edx),epsilon);
    result = [recResult1(:,1:length(recResult1)-1) recResult2(:,1:length(recResult2))];
else
    result = [Points(:,1) Points(:,edx)];
end
% If max distance is greater than epsilon, recursively simplify
    function d = penDistance(Pp, P1, P2)
        % find the distance between a Point Pp and a line segment between P1, P2.
        d = abs((P2(2,1)-P1(2,1))*Pp(1,1) - (P2(1,1)-P1(1,1))*Pp(2,1) + P2(1,1)*P1(2,1) - P2(2,1)*P1(1,1)) ...
            / sqrt((P2(2,1)-P1(2,1))^2 + (P2(1,1)-P1(1,1))^2);
    end
end