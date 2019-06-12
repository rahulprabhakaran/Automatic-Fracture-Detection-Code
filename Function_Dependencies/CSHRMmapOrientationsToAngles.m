function angles = CSHRMmapOrientationsToAngles( orientations,shearLevel )
% CSHRMmapOrientationsToAngles Map shear parameters to their associated angles.
% 
% Usage:
% 
%  angles = CSHRMmapOrientationsToAngles(orientations,shearLevel)
% 
% See also: CSHRMgetEdges, CSHRMgetRidges

    [~,coneh,ks] = CSHRMgetConeOris(shearLevel);

    nOris = numel(orientations);
    angles = zeros(size(orientations));


    for j = 1:nOris
        ori = orientations(j);
        cone = (ori >= min(coneh(:)) && ori <= max(coneh(:)));
        delta = ori - floor(ori);
        if floor(ori) == numel(ks)
            k = (1-delta)*ks(floor(ori));
        else
            k = (1-delta)*ks(floor(ori)) + delta*ks(floor(ori)+1);
        end
        if cone == 1
            angles(j) = pi/2 + atan((k/2^(shearLevel-2)));
        else
            if ori < min(coneh(:))
                angles(j) = pi-atan((k/2^(shearLevel-2)));
            else
                angles(j) = -atan((k/2^(shearLevel-2)));
            end
        end
    end

    angles = 180.*angles./pi;
    return;



    anglesMapping = [];
    ns = 1:(2^(shearLevel-2)+1);
    anglesMapping = pi - atan((ns-1)/2^(shearLevel-2));

    ns = (2^(shearLevel-2)+1):(3*2^(shearLevel-2)+1);
    anglesMapping = [anglesMapping,pi/2 - atan((ns-2^(shearLevel-1)-1)/2^(shearLevel-2))];

    ns = (3*2^(shearLevel-2)+1):2^shearLevel;
    anglesMapping = [anglesMapping,atan(-(ns-2^(shearLevel)-1)/2^(shearLevel-2))];

    anglesMapping = 180*anglesMapping/pi;
    anglesMapping = [0,anglesMapping];


    angles = anglesMapping(orientations+1);
    
end

%  Copyright (c) 2016. Rafael Reisenhofer
%
%  Part of CoShREM Toolbox v1.1
%  Built Mon, 11/01/2016
%  This is Copyrighted Material