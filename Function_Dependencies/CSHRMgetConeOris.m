function [ conev,coneh,ks ] = CSHRMgetConeOris(level)
% CSHRMgetConeOris Compute shear-indexes for the vertical and horizontal cones of a shearlet system for a specific shear-level.
% 
% Usage:
% 
%  [ conev,coneh,ks ] = CSHRMgetConeOris(level)

    ks = zeros(1,2^level+2);
      
    newOris = 1:(2^(max(level-2,0))+1);
    ks(newOris) = newOris-newOris(1);
    conev = newOris;
    
    newOris = (conev(end)+1):(conev(end)+2^(level-1)+1);
    ks(newOris) = -newOris+newOris(floor(length(newOris)/2)+1);
    coneh =  newOris;
    
    newOris = (coneh(end)+1):(2^level+2);
    if ~isempty(newOris)
        ks(newOris) = newOris-newOris(end)-1;
        conev = [conev,newOris];
    end;
    
end

%  Copyright (c) 2016. Rafael Reisenhofer
%
%  Part of CoShREM Toolbox v1.1
%  Built Mon, 11/01/2016
%  This is Copyrighted Material