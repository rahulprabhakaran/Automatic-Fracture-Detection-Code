function R = worldfileread(worldFileName, coordinateSystemType, rasterSize)
%WORLDFILEREAD Read world file and return referencing object or matrix
%
%   R = WORLDFILEREAD(worldFileName, coordinateSystemType, rasterSize)
%   reads the world file, worldFileName, and constructs a spatial
%   referencing object, R. The type of referencing object is determined by
%   the coordinateSystemType string, which can be either 'planar'
%   (including projected map coordinate systems) or 'geographic' (for
%   latitude-longitude systems). The rasterSize input should match to the
%   size of the image corresponding to the world file.
%
%   REFMAT = WORLDFILEREAD(worldFileName) reads the world file, 
%   worldFileName, and constructs a 3-by-2 referencing matrix, REFMAT.
%
%   Example 1
%   ---------
%   % Read ortho image referenced to a projected coordinate system
%   % (Massachusetts State Plane Mainland)
%   filename = 'concord_ortho_w.tif';
%   [X, cmap] = imread(filename);
%   worldFileName = getworldfilename(filename);
%   R = worldfileread(worldFileName, 'planar', size(X))
%
%   Example 2
%   ---------
%   % Read image referenced to a geographic coordinate system
%   filename = 'boston_ovr.jpg';
%   RGB = imread(filename);
%   worldFileName = getworldfilename(filename);
%   R = worldfileread(worldFileName, 'geographic', size(RGB))
%
%   See also GEOREFCELLS, GETWORLDFILENAME, MAPREFCELLS, WORLDFILEWRITE

% Copyright 1996-2016 The MathWorks, Inc.

narginchk(1, 3)
if nargin == 2
    error(message('map:validate:expected1Or3Inputs', 'WORLDFILEREAD'))
end

if nargin == 3
    % Validate coordinateSystemType, but let raster referencing functions
    % validate raster size.
    coordinateSystemType = validatestring( ...
        coordinateSystemType, {'geographic', 'planar'}, ...
        'WORLDFILEREAD', 'coordinateSystemType', 2);
end

% Check that worldFileName is a file and that it can be opened.
worldFileName = internal.map.checkfilename(worldFileName, {}, mfilename, 1, false);

% Open the input worldFileName.
fid = fopen(worldFileName);
clean = onCleanup(@() fclose(fid));

% Import the 6 world file elements into a 2-by-3 cell matrix of strings.
datastrs = scanWorldFile(fid);

% Convert the world file matrix to a raster reference object or matrix.
if nargin == 3
    if strcmp(coordinateSystemType,'geographic')
        R = worldfileToGeographicCellsReference(datastrs, rasterSize);
    else
        R = worldfileToMapCellsReference(datastrs, rasterSize);
    end
else
    R = worldfileToReferencingMatrix(datastrs);
end
end


function datastrs = scanWorldFile(fileID)
% Scan the world file, returning a 2-by-3 cell matrix of numeric data strings.
    t = textscan(fileID,'%s',6);
    map.internal.assert(numel(t{1}) == 6, 'map:fileio:expectedSixNumbers');
    datastrs = reshape(t{1},[2 3]);
end


function R = worldfileToGeographicCellsReference(datastrs, rasterSize)
% Construct a scalar map.rasterref.GeographicCellsReference object R.

    W21 = str2double(datastrs{2,1});
    W12 = str2double(datastrs{1,2});
    if W21 ~= 0 || W12 ~= 0
        error(message('map:validate:expectedRectilinearGeoRaster','WORLDFILEREAD'))
    end
    
    [deltaLatNumerator, deltaLatDenominator] = str2rat(datastrs{2,2});
    [deltaLonNumerator, deltaLonDenominator] = str2rat(datastrs{1,1});

    [num, den] = str2rat(datastrs{1,3});
    firstCornerLon = (2*num/den - deltaLonNumerator/deltaLonDenominator)/2;

    [num, den] = str2rat(datastrs{2,3});
    firstCornerLat = (2*num/den - deltaLatNumerator/deltaLatDenominator)/2;

    R = map.rasterref.GeographicCellsReference(rasterSize, ...
        firstCornerLat, firstCornerLon, ...
        deltaLatNumerator, deltaLatDenominator, ...
        deltaLonNumerator, deltaLonDenominator);
end


function R = worldfileToMapCellsReference(datastrs, rasterSize)
% Construct a scalar map.rasterref.MapCellsReference object R.

    jacobianNumerator = zeros(2,2);
    jacobianDenominator = ones(2,2);
    for k = 1:4
        [jacobianNumerator(k), jacobianDenominator(k)] = str2rat(datastrs{k});
    end
    J = jacobianNumerator ./ jacobianDenominator;
    
    [num, den] = str2rat(datastrs{1,3});
    firstCornerX = ((2*num/den) - J(1,1) - J(1,2))/2;
    
    [num, den] = str2rat(datastrs{2,3});
    firstCornerY = ((2*num/den) - J(2,1) - J(2,2))/2;

    R = map.rasterref.MapCellsReference(rasterSize, firstCornerX, ...
       firstCornerY, jacobianNumerator, jacobianDenominator);
end


function refmat = worldfileToReferencingMatrix(datastrs)
    W = zeros(2,3);
    for k = 1:numel(W)
        W(k) = str2double(datastrs{k});
    end
    refmat = worldFileMatrixToRefmat(W);
end
