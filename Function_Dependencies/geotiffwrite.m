function geotiffwrite(varargin)
%GEOTIFFWRITE Write GeoTIFF file
%
%   GEOTIFFWRITE(FILENAME, A, R) writes a georeferenced image or data grid,
%   A, spatially referenced by R, into an output file, FILENAME.
%
%   FILENAME is a character string that specifies the output file name and
%   location. If your FILENAME includes an extension, it must be '.tif' or
%   '.TIF'. The output file is a tiled GeoTIFF file if the input, A, is at
%   least 160-by-160 in size. Otherwise, it is organized as rows-per-strip.
%
%   A is an M-by-N array (grayscale image or data grid) or an M-by-N-by-P
%   array (color or hyperspectral image, or data grid). The coordinates of
%   A are geographic and in the 'WGS 84' coordinate system, unless you
%   specify 'GeoKeyDirectoryTag' or 'CoordRefSysCode' and indicate a
%   different coordinate system.
%
%   In this first syntax, R is always a geographic raster reference object,
%   a referencing matrix, or a referencing vector. (If you are working with
%   image coordinates in a projected coordinate system and R is a
%   map raster reference object or a referencing matrix, specify
%   'GeoKeyDirectoryTag' or 'CoordRefSysCode' accordingly. See the
%   Name-Value Pair Arguments section for more information.)
%
%   GEOTIFFWRITE(FILENAME, X, CMAP, R) writes the indexed image in X and
%   its associated colormap, CMAP, to FILENAME. X is spatially referenced
%   by R.
%
%   GEOTIFFWRITE(...,NAME1, VALUE1, NAME2, VALUE2,...) writes a
%   georeferenced image or data grid with additional options that control
%   various characteristics of the output file specified by one or more
%   NAME, VALUE pair arguments. NAME is the argument name and VALUE is the
%   corresponding value. NAME must appear inside single quotes ('') and is
%   case insensitive. You can specify several name-value pair arguments in
%   any order.
%
%   Class Support
%   -------------
%   The input array, A, can be any numeric class or logical.
%
%   The indexed image, X, must be class uint8 or uint16. The colormap
%   array, CMAP, must be class double.
%
%   Name-Value Pair Arguments
%   -------------------------
%   'CoordRefSysCode'        Scalar, positive, integer-valued number that 
%                            specifies the coordinate reference system code
%                            for the coordinates of the data.
%
%                            You can specify coordinates in either a
%                            geographic or a projected coordinate system,
%                            and you can use a string, such as,
%                            'EPSG:4326'. If you specify the coordinate
%                            system with a string, include the 'EPSG:'
%                            prefix. See the GeoTIFF specification or the
%                            EPSG data files (pcs.csv and gcs.csv) for the
%                            code numbers.
%
%                            If you specify both the 'GeoKeyDirectoryTag'
%                            and the 'CoordRefSysCode', the coordinate
%                            system code found in the 'CoordRefSysCode'
%                            takes precedence over the coordinate system
%                            key found in the 'GeoKeyDirectoryTag'. If one
%                            value specifies a geographic coordinate system
%                            and the other value specifies a projected
%                            coordinate system, you receive an error.
%
%                            If you do not specify a value for this
%                            argument, the default value is 4326,
%                            indicating that the coordinates are geographic
%                            and in the 'WGS 84' geographic coordinate
%                            system.
%
%   'GeoKeyDirectoryTag'     Structure that specifies the GeoTIFF
%                            coordinate reference system and
%                            meta-information. The structure contains field
%                            names that match the GeoKey names found in the
%                            GeoTIFF specification. The field names are
%                            case insensitive. The structure can be
%                            obtained from the GeoTIFF information
%                            structure, returned by the function 
%                            <a href="matlab:doc geotiffinfo">geotiffinfo</a>, in the field, 
%                            GeoTIFFTags.GeoKeyDirectoryTag.
%                            If you set certain fields of the
%                            'GeoKeyDirectoryTag' to inconsistent settings,
%                            you receive an error message. See the
%                            <a href="matlab:doc geotiffwrite">geotiffwrite</a> reference page 
%                            for more details about error messages.
%
%   'RPCCoefficientTag'      Scalar map.geotiff.RPCCoefficientTag that
%                            specifies an optional RPCCoefficientTag.
%
%   'TiffTags'               Structure that specifies values for the TIFF
%                            tags in the output file. The field names of
%                            the structure match the TIFF tag names
%                            supported by the Tiff class.  The field names
%                            are case insensitive.  See the
%                            <a href="matlab:doc geotiffwrite">geotiffwrite</a> reference page 
%                            for details about 'TiffTags'.
%
%   References
%   ----------
%   The 'CoordRefSysCode' value for geographic coordinate systems may be
%   found in the GeoTIFF specification at:
%   <a href="matlab:web('http://geotiff.maptools.org/spec/geotiff6.html#6.3.2.1')
%   ">http://geotiff.maptools.org/spec/geotiff6.html#6.3.2.1</a>
%
%   The 'CoordRefSysCode' value for projected coordinate systems may be
%   found in the GeoTIFF specification at:
%   <a href="matlab:web('http://geotiff.maptools.org/spec/geotiff6.html#6.3.3.1')
%   ">http://geotiff.maptools.org/spec/geotiff6.html#6.3.3.1</a>
% 
%   The GeoKey field names for the 'GeoKeyDirectoryTag' may be found in the
%   GeoTIFF specification in the 6.2 Key ID Summary section at:
%   <a href="matlab:web('http://geotiff.maptools.org/spec/geotiff6.html#6.2')
%   ">http://geotiff.maptools.org/spec/geotiff6.html#6.2</a>
%
%   The 'CoordRefSysCode' values may also be obtained from the EPSG data
%   files (pcs.csv and gcs.csv) in the directory
%   <matlabroot>/toolbox/map/mapproj/projdata/epsg_csv
%
%   Example 1
%   ---------
%   % Write an image with geographic coordinates from a JPEG file to a
%   % GeoTIFF file.
%   basename = 'boston_ovr';
%   imagefile = [basename '.jpg'];
%   RGB = imread(imagefile);
%   worldfile = getworldfilename(imagefile);
%   R = worldfileread(worldfile, 'geographic', size(RGB));
%   filename = [basename '.tif'];
%   geotiffwrite(filename, RGB, R)
%   figure
%   usamap(RGB, R)
%   geoshow(filename)
%
%   Example 2
%   ---------
%   % Write a WMS image to a GeoTIFF file.
%   nasa = wmsfind('nasa', 'SearchField', 'serverurl');
%   layerName = 'bluemarbleng';
%   layer = nasa.refine(layerName,  'SearchField', 'layername', ...
%        'MatchType', 'exact');
%   [A, R] = wmsread(layer(1));
%   filename = [layerName '.tif'];
%   geotiffwrite(filename, A, R)
%   figure
%   worldmap world
%   geoshow(filename)
%
%   Example 3
%   ---------
%   % Write the Concord orthophotos to a single GeoTIFF file.
%   % Read the two adjacent orthophotos and combine them.
%   X_west = imread('concord_ortho_w.tif');
%   X_east = imread('concord_ortho_e.tif');
%   X = [X_west X_east];
%
%   % Construct referencing objects for the orthophotos and for their
%   % combination.
%   R_west = worldfileread('concord_ortho_w.tfw', 'planar', size(X_west));
%   R_east = worldfileread('concord_ortho_e.tfw', 'planar', size(X_east));
%   R = R_west;
%   R.XWorldLimits = [R_west.XWorldLimits(1) R_east.XWorldLimits(2)];
%   R.RasterSize = size(X);
%
%   % Write the combined image to a GeoTIFF file. Use the code number, 
%   % 26986, indicating the PCS_NAD83_Massachusetts Projected Coordinate 
%   % System.
%   coordRefSysCode = 26986;
%   filename = 'concord_ortho.tif';
%   geotiffwrite(filename, X, R, 'CoordRefSysCode', coordRefSysCode);
%   figure
%   mapshow(filename)
%
%   Example 4
%   ---------
%   % Write the first 1024 columns and last 1024 rows of a GeoTIFF file 
%   % to a new GeoTIFF file.
%   [A, R] = geotiffread('boston.tif');
%   row = [size(A,1)-1024+1 size(A,1)];
%   col = [1 1024];
%   subImage = A(row(1):row(2), col(1):col(2), :);
%   xi = col + [-.5 .5];
%   yi = row + [-.5 .5];
%   [xlimits, ylimits] = intrinsicToWorld(R, xi, yi);
%   subR = R;
%   subR.RasterSize = size(subImage);
%   subR.XWorldLimits = sort(xlimits);
%   subR.YWorldLimits = sort(ylimits);
%   info = geotiffinfo('boston.tif');
%   filename = 'boston_subimage.tif';
%   geotiffwrite(filename, subImage, subR,  ...
%      'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
%   figure
%   mapshow(filename);
%
%   Example 5
%   ---------
%   % Write the Mount Washington SDTS DEM terrain data to GeoTIFF.
%   % The data are referenced to Universal Transverse Mercator (UTM), 
%   % Zone 19, in the North American Datum of 1927. This corresponds
%   % to the GeoTIFF PCS_NAD27_UTM_zone_19N code number 26719. Set the 
%   % raster interpretation to 'postings' because it's a USGS DEM. This
%   % corresponds to the GeoTIFF raster type PixelIsPoint.
%   [Z, refmat] = sdtsdemread('9129CATD.ddf');
%   R = refmatToMapRasterReference(refmat, size(Z),'postings');
%   key.GTModelTypeGeoKey  = 1;  % Projected Coordinate System (PCS)
%   key.GTRasterTypeGeoKey = 2;  % PixelIsPoint
%   key.ProjectedCSTypeGeoKey = 26719;
%   filename = '9129.tif';
%   geotiffwrite(filename, Z, R, 'GeoKeyDirectoryTag', key);
%
%   % Plot the outline of the state of New Hampshire in UTM.
%   S = shaperead('usastatelo', 'UseGeoCoords', true, 'Selector',...
%      {@(name) any(strcmp(name,{'New Hampshire'})), 'Name'});
%   proj = geotiffinfo(filename);
%   [x, y] = projfwd(proj, [S.Lat], [S.Lon]);
%   figure
%   mapshow(x,y)
%
%   % Display the GeoTIFF DEM file.
%   hold on
%   h = mapshow(filename, 'DisplayType', 'surface');
%   demcmap(get(h,'ZData'))
%   
%   See also GEOTIFFINFO, GEOTIFFREAD, IMREAD, IMWRITE, map.geotiff.RPCCoefficientTag, Tiff

% Copyright 2010-2017 The MathWorks, Inc.

% Verify the number of inputs.
narginchk(3,inf);

% Parse the inputs.
[filename, A, cmap, R, Params] = parseInputs(varargin{:});

% Validate the inputs.
[filename, A, cmap, R, Params] = validateInputs(filename, A, cmap, R, Params);

% Construct TIFF tags structure based on input image characteristics and
% the tag values in the TiffTags parameter.
TiffTags = constructTiffTags(A, cmap, Params.TiffTags);

% Construct GeoTIFF tags structure based on the spatial referencing
% information from R and the coordinate reference information in
% the GeoKeyDirectoryTag parameter.
GeoTiffTags = constructGeoTiffTags(R, Params.GeoKeyDirectoryTag);

% Write the GeoTIFF file.
writeGeoTiffFile(filename, A, TiffTags, GeoTiffTags, Params.RPCCoefficientTag);

%--------------------------------------------------------------------------

function [filename, A, cmap, R, Params] = parseInputs(filename, A, varargin)
% Parse the input arguments. Prior to calling this function, VARARGIN must
% be validated to contain at least one element (or the calling workspace
% varargin must be validated to contain at least three elements.) The
% values or expected class of the arguments are not validated in this
% function.
% 
% FILENAME and A are the first two elements of VARARGIN (from the calling
% workspace) and can be assigned directly using input variable names.
%
% PARAMS is a structure containing field names that match the parameter
% names of the interface:
%  'TiffTags', 'CoordRefSysCode', 'GeoKeyDirectoryTag', 'RPCCoefficientTag'
% The fields are empty if the parameter is not supplied in VARARGIN.

% Get number of data parameters.
numDataArgs = internal.map.getNumberOfDataArgs(varargin{:});

% Obtain required data parameters.
switch numDataArgs
    % (FILENAME, A, 'PARAM1', ...)
    case 0
        error(message('map:geotiff:tooFewDataArgs', 'GEOTIFFWRITE', 'A', 'R'))
        
    case 1
        % (FILENAME, A, R, 'PARAM1', ...)
        R = varargin{1};
        cmap = [];
        varargin(1) = [];
              
    case 2
        % (FILENAME, X, CMAP, R, 'PARAM1' ...)       
        cmap = varargin{1};
        R = varargin{2};
        varargin(1:2) = [];
        
    otherwise
        % Include filename and A in the count.
        error(message('map:geotiff:tooManyDataArgs', 'GEOTIFFWRITE', 'A', 'CMAP', 'R', numDataArgs + 2))
end

% Parse the optional parameters and return a structure.
Params = parseParams(varargin);

%--------------------------------------------------------------------------

function Params = parseParams(inputs)
% Parse the optional input parameters and return them in the structure
% PARAMS which contain field names that match the parameter names of the
% interface:
%  'TiffTags', 'CoordRefSysCode', 'GeoKeyDirectoryTag', 'RPCCoefficientTag'
% The fields are empty [] if the parameter is not supplied in the cell 
% array INPUTS.

map.internal.assert(mod(numel(inputs), 2) == 0, ...
    'map:geotiff:parameterCountNotEven');

% Create output structure with paramNames as the field names of the
% structure. The field values are all [].
paramNames = {'TiffTags', 'CoordRefSysCode', 'GeoKeyDirectoryTag', 'RPCCoefficientTag'};
Params = cell2struct(cell(size(paramNames)), paramNames, 2);

% Copy the 'PARAM1', VALUE1, ... parameters to the Params structure. Use
% validatestring to match the field names with the parameter names.
if ~isempty(inputs)
    names  = inputs(1:2:end);
    values = inputs(2:2:end);
    
    for k=1:numel(names)
        name = validatestring(names{k}, paramNames);
        Params.(name) = values{k};
    end
end

%--------------------------------------------------------------------------

function [filename, A, cmap, R, Params] = validateInputs( ...
    filename, A, cmap, R, Params)
% Validate the inputs.

% Validate filename.
filename = validateFilename(filename);

% Validate image (or data grid) parameters.
[A, cmap] = validateImage(A, cmap);

% Validate Params.
Params = validateParams(Params);

% Obtain the coordinate system type: 'projected' or 'geographic'
type = getGeoKeyCRSType(Params.GeoKeyDirectoryTag);

% Validate R and if required, convert it to a referencing object.
hasColorMap = ~isempty(cmap);
R = validateR(R, size(A), hasColorMap, type);

% Validate GeoKeyDirectoryTag raster type code ('PixelIsArea' or
% 'PixelIsPoint') and update the GeoKeyDirectoryTag based on the raster
% type, if needed.
Params.GeoKeyDirectoryTag = validateRasterType(R, Params.GeoKeyDirectoryTag);

% Validate image shape (rows/columns start from position).
[A, R] = validateImageShape(A, R);

% Validate the photometric interpretation value, if set, with the image.
validatePhotometricWithImage(Params.TiffTags, A);

%--------------------------------------------------------------------------

function filename = validateFilename(filename)
% Validate the filename input.

validateattributes(filename, {'char'}, {'vector','nonempty'}, ...
    mfilename, 'FILENAME', 1);

% Add .tif extension to filename if necessary.
tifExt = '.tif';
[pathname, basename, ext] = fileparts(filename);
map.internal.assert(strcmpi(ext, tifExt) || isempty(ext), ...
   'map:geotiff:invalidExtension', upper( mfilename ), filename, tifExt);
filename = fullfile(pathname,[basename, tifExt]);

%--------------------------------------------------------------------------

function [A, cmap] = validateImage(A, cmap)
% Validate the input image A (or data grid) and its associated colormap
% (which is empty if the image is not indexed).

if isempty(cmap)
    % M-by-N, M-by-N-by-P (data grid or image)
    validateattributes(A, {'numeric', 'logical'}, {'nonempty'}, ...
        mfilename, 'A', 2);
else
    % M-by-N indexed image
    validateattributes(A, {'uint8','uint16'}, {'nonempty'}, ...
        mfilename, 'X', 2);
    internal.map.checkcmap(cmap, mfilename, 'CMAP', 3)
    
    % Pad the colormap if needed.
    classType = class(A);
    maxNumColors = double(intmax(classType)) + 1;
    numRows = size(cmap,1);
    map.internal.assert(numRows <= maxNumColors, ...
       'map:geotiff:tooManyColormapValues', maxNumColors, numRows);
    cmap(end+1:maxNumColors,:) = 0;
end

%--------------------------------------------------------------------------

function Params = validateParams(S)
% Validate the optional parameter/value pair structure. A new structure is 
% returned containing only the field names:
%    'TiffTags', 'GeoKeyDirectoryTag', and 'RPCCoefficientTag'
%
% The 'CoordRefSysCode' field of S is merged into the 'GeoKeyDirectoryTag'
% field.

% Validate TiffTags field.
Params.TiffTags = validateTiffTags(S.TiffTags);

% Validate CoordRefSysCode field.
coordRefSysCode = validateCoordRefSysCode(S.CoordRefSysCode);

% Validate GeoKeyDirectoryTag field.
GeoKeyDirectoryTag = validateGeoKeyDirectoryTag(S.GeoKeyDirectoryTag);

% Merge CoordRefSysCode into the GeoKeyDirectoryTag.
Params.GeoKeyDirectoryTag = mergeCoordRefSysCode( ...
    GeoKeyDirectoryTag, coordRefSysCode);

% Validate RPCCoefficientTag, if present.
if ~isempty(S.RPCCoefficientTag)
    validateattributes(S.RPCCoefficientTag, {'map.geotiff.RPCCoefficientTag'}, ...
        {'scalar'}, mfilename, 'RPCCoefficientTag')
end
Params.RPCCoefficientTag = S.RPCCoefficientTag;

%--------------------------------------------------------------------------

function TiffTags = validateTiffTags(TiffTags)
% Validate the 'TiffTags' parameter value.

if isempty(TiffTags)
    TiffTags = struct();
else
    % Validate TiffTags as a scalar structure.
    validateattributes( ...
        TiffTags, {'struct'}, {'scalar'}, mfilename, 'TiffTags')
    
    % Validate the field names of TiffTags.
    TiffTags = validateTiffTagNames(TiffTags);
    
    % Validate Compression field and change string value to numeric value.
    if isfield(TiffTags, 'Compression')
        TiffTags.Compression = validateCompression(TiffTags.Compression);
    end
    
    % Validate Photometric field and change string value to numeric value.
    if isfield(TiffTags, 'Photometric')
        TiffTags.Photometric = validatePhotometric(TiffTags.Photometric);
    end
    
    % Validate RowsPerStrip with TileLength or TileWidth.
    if isfield(TiffTags, 'RowsPerStrip') && ...
            (isfield(TiffTags, 'TileLength') || ...
             isfield(TiffTags, 'TileWidth'))
         error(message('map:geotiff:invalidImageOrganization', 'TiffTags', 'RowsPerStrip', 'TileLength', 'TileWidth'));
    end
end

%--------------------------------------------------------------------------

function TiffTags = validateTiffTagNames(TiffTags)
% Validate the field names of TiffTags to be valid Tiff tag names excluding
% the reserved names. The reserved names are either:
%    1) set by the function,
%    2) GeoTIFF tag names, or
%    3) read-only tag names.

% Obtain all the tag names from the Tiff class.
tiffTagNames = Tiff.getTagNames;

% Replace ASCII with Ascii in order to be consistent in any generated error
% messages (GEOTIFFINFO uses GeoAsciiParamsTag).
tiffTagNames = strrep(tiffTagNames,'GeoASCIIParamsTag','GeoAsciiParamsTag');

% Remove 'Photometric' and replace with 'PhotometricInterpretation'
% validatestring will still allow 'Photometric' to be used with this
% substitution. Replace the field after using validatestring.
tiffTagNames = strrep(tiffTagNames,'Photometric','PhotometricInterpretation');

% Assign the names for the reserved tags.
reservedTagNames = { ...
    'BitsPerSample', 'SampleFormat', 'SamplesPerPixel', ...
    'StripByteCounts', 'StripOffsets', 'SubFileType', 'SubIFD', ...
    'TileByteCounts', 'TileOffsets', ...
    'ColorMap', 'ImageLength', 'ImageWidth', ...
    'GeoAsciiParamsTag', 'GeoDoubleParamsTag', 'GeoKeyDirectoryTag', ...
    'ModelPixelScaleTag', 'ModelTiepointTag', 'ModelTransformationTag'};

% Validate the field names, allowing only the permitted names to be set.
permittedNames = setdiff(tiffTagNames, reservedTagNames);
tagNames = fieldnames(TiffTags);
for k = 1:numel(tagNames)
    inputName = tagNames{k};
    tagName = validatestring( ...
        inputName, permittedNames, mfilename, ['TiffTags.' inputName]);
    if ~isfield(TiffTags, tagName)
        % tagName and inputName may not match even though a match is found
        % since the string may be a partial match or not match case.
        % Switch the case and/or change partial/lower case to mixed
        % case to match the field names of the Tiff class.
        TiffTags.(tagName) = TiffTags.(inputName);
        TiffTags = rmfield(TiffTags, inputName);
    end
end

% Replace 'PhotometricInterpretation' with 'Photometric' (for Tiff class)
photometricInterpretation = 'PhotometricInterpretation';
if isfield(TiffTags, photometricInterpretation)
    TiffTags.Photometric = TiffTags.(photometricInterpretation);
    TiffTags = rmfield(TiffTags, photometricInterpretation);
end

%--------------------------------------------------------------------------

function compression = validateCompression(compression)
% Validate the input COMPRESSION. If it is a string, convert the expected
% value to the value expected by the Tiff class.

% Modify 'uncompressed' to 'none' for compatibility with imfinfo.
if strcmpi(compression, 'uncompressed')
    compression = 'none';
end

permittedValues = {'LZW', 'PackBits', 'Deflate', 'None'};
compression = tiffNameToValue(compression, permittedValues, 'Compression');

%--------------------------------------------------------------------------

function photometric = validatePhotometric(photometric)
% Validate the input PHOTOMETRIC. If it is as a string, convert the
% expected string value to the value expected by the Tiff class.

permittedValues = {'MinIsBlack', 'RGB', 'Palette', 'Separated'};
photometric = tiffNameToValue(photometric, permittedValues, 'Photometric');

%--------------------------------------------------------------------------

function value = tiffNameToValue(value, permittedValues, fieldValue)
% Validate the input, VALUE, to be an element in the cell array,
% permittedValues or a numeric value that is defined by the property,
% Tiff.(fieldValue). If VALUE is a string, convert the valid string to the
% numeric value that is expected by the Tiff class (Tiff.(fieldValue))

% Obtain the values assigned by the Tiff class for the permittedValues.
values = zeros(numel(permittedValues), 1);
for k=1:numel(values)
    values(k) = Tiff.(fieldValue).(permittedValues{k});
end

% Obtain the numerical value.
valueName = ['TiffTags.' fieldValue];
if ischar(value)
    strvalue = validatestring(value, permittedValues, mfilename, valueName);
    index = strcmp(permittedValues, strvalue);
    value = values(index);
else
    validateattributes(value, ...
        {'numeric'}, {'scalar', 'finite'},  mfilename, valueName) 
    index = value == values;
    if ~any(index)
        cvalues = sprintf('%d, ', values);
        cvalues(end-1:end) = '. ';
        error(message('map:geotiff:invalidTiffTagValue', valueName, cvalues, value))
    end
end
      
%--------------------------------------------------------------------------

function GeoKeyTag = validateGeoKeyDirectoryTag(GeoKeyTag)
% Validate the 'GeoKeyDirectoryTag' parameter value.

if isempty(GeoKeyTag)
    GeoKeyTag = struct();
else
    % Validate as a scalar structure.
    validateattributes(GeoKeyTag, ...
        {'struct'}, {'scalar'}, mfilename, 'GeoKeyDirectoryTag')
    
    % Validate field names.
    GeoKeyTag = validateGeoKeyDirectoryTagNames(GeoKeyTag);
    
    % Validate data type of field values.
    validateGeoKeyDirectoryTagValues(GeoKeyTag);
    
    % Validate the coordinate reference system if defined in the
    % GeoKeyDirectoryTag value.
    GeoKeyTag = validateGeoKeyCRS(GeoKeyTag);
end

%--------------------------------------------------------------------------

function GeoKeyTag = validateGeoKeyDirectoryTagNames(GeoKeyTag)
% Validate the field names of the GeoKeyDirectoryTag.

% Construct the key ID to name map.
idToNameMap = constructGeoKeyDirectoryMap;
permittedNames = idToNameMap.values;

% Remove 'UserDefined'
permittedNames = strrep(permittedNames, 'UserDefined', '');

keyNames = fieldnames(GeoKeyTag);
for k = 1:numel(keyNames)
    inputName = keyNames{k};
    keyName = validatestring(inputName, permittedNames, mfilename,  ...
        ['GeoKeyDirectoryTag.' inputName]);
    if ~isfield(GeoKeyTag, keyName)
        % Switch case or partial name.
       GeoKeyTag.(keyName) = GeoKeyTag.(inputName);
       GeoKeyTag = rmfield(GeoKeyTag, inputName);
    end
end

%--------------------------------------------------------------------------

function validateGeoKeyDirectoryTagValues(GeoKeyTag)
% Validate the values in the 'GeoKeyDirectoryTag' structure, GeoKeyTag. 

% Any *Citation* key name must be a char and not empty.
names = fieldnames(GeoKeyTag);
cellIndex = regexpi(names, regexptranslate('wildcard', 'Citation'));
stringFields = names(~cellfun(@isempty, cellIndex));
numericFields = setdiff(names, stringFields);
for k=1:numel(stringFields)
    field = stringFields{k};
    validateattributes(GeoKeyTag.(field), ...
        {'char'}, {'nonempty'}, mfilename, ['GeoKeyDirectoryTag.' field])
end

% Validate other fields as numeric, non-empty, and finite.
for k=1:numel(numericFields)
    field = numericFields{k};
    validateattributes(GeoKeyTag.(field), {'numeric'},  ...
        {'nonempty', 'finite'}, mfilename, ['GeoKeyDirectoryTag.' field])
end

%--------------------------------------------------------------------------

function GeoKeyTag = validateGeoKeyCRS(GeoKeyTag)
% Validate consistency of coordinate reference system information in the
% 'GeoKeyDirectoryTag' structure, GeoKeyTag.

modelField = 'GTModelTypeGeoKey';
projField  = 'ProjectedCSTypeGeoKey';
geoField   = 'GeographicTypeGeoKey';

tfModelField = isfield(GeoKeyTag, modelField);
tfProjField  = isfield(GeoKeyTag, projField);
tfGeoField   = isfield(GeoKeyTag, geoField);

if tfProjField && tfGeoField
    % User has specified both geographic and projected systems, which is
    % valid. However, additional information must be specified.
    % Expect the GTModelTypeGeoKey to defined and to indicate projected
    % coordinates. (The coordinate reference system code is validated in the
    % next if block).
    map.internal.assert(tfModelField && GeoKeyTag.(modelField) == 1, ...
       'map:geotiff:geoAndMapProjections', projField, geoField, modelField);
end

if tfProjField
    % projected
    key = GeoKeyTag.(projField);
    geoKeyType = internal.map.getCoordRefSysCodeType(key);
    % The geoKeyType must be 'projected' or user-defined.
    map.internal.assert(any(strcmp(geoKeyType, {'projected','user-defined'})), ...
       'map:geotiff:invalidProjectedCSTypeGeoKey', [ 'GeoKeyDirectoryTag.', projField ])
    
    if tfModelField
        % GTModelTypeGeoKey = 1 indicates projected coordinates
        map.internal.assert(GeoKeyTag.(modelField) == 1, ...
           'map:geotiff:invalidProjectedModelType', [ 'GeoKeyDirectoryTag.', modelField ])
    else
        % The GTModelTypeGeoKey is not required to be set by the user but
        % is required to be set in the file. Since it is known now, 
        % set the value to indicate projected coordinates.
        GeoKeyTag.(modelField) = 1;
    end
    
elseif tfGeoField
    % geographic
    key = GeoKeyTag.(geoField);
    geoKeyType = internal.map.getCoordRefSysCodeType(key);
    % The geoKeyType must be 'geographic' or user-defined.
    map.internal.assert(any(strcmp(geoKeyType, {'geographic', 'user-defined'})), ...
      'map:geotiff:invalidGeographicTypeGeoKey', [ 'GeoKeyDirectoryTag.', geoField ])
    
    if tfModelField
        % GTModelTypeGeoKey = 2 indicates geographic coordinates
        map.internal.assert(GeoKeyTag.(modelField) == 2, ...
           'map:geotiff:invalidGeographicModelType', [ 'GeoKeyDirectoryTag.', modelField ])
    else
        % The GTModelTypeGeoKey is not required to be set by the user but
        % is required to be set in the file. Since it is known now, 
        % set the value to indicate geographic coordinates.
        GeoKeyTag.(modelField) = 2;
    end
end

%--------------------------------------------------------------------------

function coordRefSysCode = validateCoordRefSysCode(coordRefSysCode)
% Validate the 'CoordRefSysCode' parameter value. Return a string valued
% 'EPSG:' code to its equivalent double value.

if ~isempty(coordRefSysCode)   
    strName = 'CoordRefSysCode';
    if ischar(coordRefSysCode)
        value = coordRefSysCode;
        epsg = 'EPSG:';
        map.internal.assert(strncmpi(value, epsg, numel(epsg)), ...
           'map:geotiff:notEPSG', strName, epsg)
        coordRefSysCode = str2double(value(numel(epsg)+1:end));
    end
    
    % Validate class type.
    validateattributes(coordRefSysCode, ...
        {'numeric'}, {'scalar', 'finite'},  mfilename, strName)
  
    % Validate coordinate type: 'projected' or 'geographic'.
    coordRefSysCodeType = internal.map.getCoordRefSysCodeType(coordRefSysCode);
    map.internal.assert(strcmp(coordRefSysCodeType, 'projected') ...
        || strcmp(coordRefSysCodeType, 'geographic'), ...
        'map:geotiff:invalidCoordRefSysCodeType', 'CoordRefSysCode');
end

%--------------------------------------------------------------------------

function GeoKeyTag = mergeCoordRefSysCode(GeoKeyTag, coordRefSysCode)
% Merge the coordRefSysCode into the GeoKeyTag structure.

geoField  = 'GeographicTypeGeoKey';
projField = 'ProjectedCSTypeGeoKey';
tfProjField = isfield(GeoKeyTag, projField);
tfGeoField  = isfield(GeoKeyTag, geoField);

if ~tfProjField && ~tfGeoField && isempty(coordRefSysCode)
    % Set default value for coordinate reference system.
    coordRefSysCode = 4326;
    GeoKeyTag.(geoField) = coordRefSysCode;
end

% 'GeoKeyDirectoryTag' coordinate type.
geoKeyType = getGeoKeyCRSType(GeoKeyTag);

% 'CoordRefSysCode'
if ~isempty(coordRefSysCode)
    coordRefSysCodeType = internal.map.getCoordRefSysCodeType(coordRefSysCode);
    if ~isempty(geoKeyType)
        % Validate that the codes represent the same type of
        % coordinate system.
        map.internal.assert(strcmp(coordRefSysCodeType, geoKeyType), ...
           'map:geotiff:multipleCoordinateReferenceSystemTypes', 'CoordRefSysCode', geoKeyType);
    end
    % Merge the CoordRefSysCode into the GeoKeyTag structure.
    GeoKeyTag = mergeCoordRefSysCodeType( ...
        GeoKeyTag, coordRefSysCode, coordRefSysCodeType);
end
   
%--------------------------------------------------------------------------

function geoKeyType = getGeoKeyCRSType(GeoKeyTag)
% Return the coordinate system type of the input GeoKeyTag structure. The
% output string, geoKeyType will either be:
% 'geographic', 'projected', or '' if undetermined.

projField  = 'ProjectedCSTypeGeoKey';
geoField   = 'GeographicTypeGeoKey';

tfProjField  = isfield(GeoKeyTag, projField);
tfGeoField   = isfield(GeoKeyTag, geoField);
    
if tfProjField
    % projected
    key = GeoKeyTag.(projField);
    geoKeyType = internal.map.getCoordRefSysCodeType(key);
    if strcmp(geoKeyType, 'user-defined')
        geoKeyType = 'projected';
    end
    
elseif tfGeoField
    % geographic
    key = GeoKeyTag.(geoField);
    geoKeyType = internal.map.getCoordRefSysCodeType(key);
    if strcmp(geoKeyType, 'user-defined')
        geoKeyType = 'geographic';
    end
    
else
    geoKeyType = '';
end

%--------------------------------------------------------------------------

function GeoKeyTag = mergeCoordRefSysCodeType(GeoKeyTag, code, type)
% Merge the 'CoordRefSysCode', code into the GeoKeyTag structure based on
% the value of the string, type, indicating either 'projected' or
% 'geographic'. Prior to calling this function, type must be set to one
% value or the other (and is not validated in this function). The GeoKeyTag
% fields are set to indicate either a projected or geographic coordinate
% system and the projection code is set.

modelField = 'GTModelTypeGeoKey';
tfModelField = isfield(GeoKeyTag, modelField);

if strcmp(type, 'projected')
    if tfModelField
        % GTModelTypeGeoKey = 1 indicates projected coordinates
        map.internal.assert(GeoKeyTag.(modelField) == 1, ...
          'map:geotiff:invalidProjectedModelType', [ 'GeoKeyDirectoryTag.', modelField ])
    else
        GeoKeyTag.(modelField) = 1;
    end
    GeoKeyTag.ProjectedCSTypeGeoKey = code;
else
    if tfModelField
        % GTModelTypeGeoKey = 2 indicates geographic coordinates
        map.internal.assert(GeoKeyTag.(modelField) == 2, ...
           'map:geotiff:invalidGeographicModelType', [ 'GeoKeyDirectoryTag.', modelField ])
    else
        GeoKeyTag.GTModelTypeGeoKey = 2;
    end
    GeoKeyTag.GeographicTypeGeoKey = code;
end

%--------------------------------------------------------------------------

function R = validateR(R, rasterSize, hasColorMap, coordRefSysType)
% Validate R and return a referencing object.

if hasColorMap
    % (FILENAME, A, CMAP, R)
    arg_pos = 4;
else
    % (FILENAME, A, R)
    arg_pos = 3;
end

if strcmp(coordRefSysType, 'geographic')
    if isobject(R) && isscalar(R) && isprop(R,'CoordinateSystemType') ...
            && strcmp(R.CoordinateSystemType,'planar')
        error(message('map:geotiff:expectedGeoRasterReference', 'R', ...
            class(R), 'CoordRefSysCode', 'GeoKeyDirectoryTag'));
    else
        angle_units = 'degrees';
        R = internal.map.convertToGeoRasterRef(R, rasterSize, angle_units, ...
            upper(mfilename), 'R', arg_pos);
    end   
elseif strcmp(coordRefSysType, 'projected')
    R = convertToMapRasterRef(R, rasterSize, upper(mfilename), 'R', arg_pos);
end
    
%--------------------------------------------------------------------------

function R = convertToMapRasterRef( ...
    R, rasterSize, func_name, var_name, arg_pos)
% Validate and convert R to a map raster reference object.

if numel(R) == 3
    error(message('map:geotiff:refvecNotAllowed','R','map raster reference'));
end

R = refmatToMapRasterReference(R, rasterSize, func_name, var_name, arg_pos);

map.internal.assert(strcmp(R.TransformationType, 'rectilinear'), ...
   'map:geotiff:skewOrRotation', func_name, func_name)

%--------------------------------------------------------------------------

function GeoKeyTag = validateRasterType(R, GeoKeyTag)
% Validate GeoKeyDirectoryTag raster type code to be consistent with the
% referencing object.

if strcmp(R.RasterInterpretation, 'cells')
    % RasterPixelIsArea
    rasterType = 1;
else
    % RasterPixelIsPoint
    rasterType = 2;
end

field = 'GTRasterTypeGeoKey';
if isfield(GeoKeyTag, field)
    map.internal.assert(rasterType == GeoKeyTag.(field), ...
       'map:geotiff:inconsistentRasterType', [ 'GeoKeyDirectoryTag.', field ],  ...
          rasterType, R.RasterInterpretation, GeoKeyTag.(field));
else
    GeoKeyTag.(field) = rasterType;
end

%--------------------------------------------------------------------------

function [A, R] = validateImageShape(A, R)
% Validate the image shape and flip the rows or columns if needed to be
% consistent with the GeoTIFF specification.

% Assuming rectilinear.
if strcmp(R.ColumnsStartFrom, 'south')
    R.ColumnsStartFrom = 'north';
    A = A(end:-1:1,:,:);
end

if strcmp(R.RowsStartFrom, 'east')
    R.RowsStartFrom = 'west';
    A = A(:,end:-1:1,:);
end

%--------------------------------------------------------------------------

function validatePhotometricWithImage(TiffTags, A)
% Validate the photometric interpretation value with the image.

if isfield(TiffTags, 'Photometric') ...
        && TiffTags.Photometric == Tiff.Photometric.RGB
    map.internal.assert(ndims(A) == 3 && size(A,3) == 3, ...
       'map:geotiff:imageSizeNotRGB', 'PhotometricInterpretation')
end
        
%--------------------------------------------------------------------------

function GeoTiffTags = constructGeoTiffTags(R, GeoKeyTag)
% Construct GeoTIFF tags from R and the GeoKeyTag structure.

% R is validated to be rectilinear. If this restriction is changed, this
% code block must be modified to write out either a ModelPixelScaleTag,
% ModelTiepointTag or ModelTransformationTag.

% ModelPixelScaleTag
pixelScale = constructModelPixelScale(R);
GeoTiffTags.ModelPixelScaleTag = pixelScale;

% ModelTiepointTag
tiepoint = constructModelTiepoint(R);
GeoTiffTags.ModelTiepointTag = tiepoint;

% GeoKeyDirectoryTag
[keys, S] = constructGeoKeyDirectory(GeoKeyTag);
GeoTiffTags.GeoKeyDirectoryTag = keys;

% GeoDoubleParamsTag
if ~isempty(S.DoubleNames)
    params = constructGeoDoubleParams(GeoKeyTag, S.DoubleNames, S.DoubleOffset);
    GeoTiffTags.GeoDoubleParamsTag = params;
end

% GeoAsciiParamsTag
if ~isempty(S.AsciiNames)
    params = constructGeoAsciiParams(GeoKeyTag, S.AsciiNames, S.AsciiOffset);
    GeoTiffTags.GeoASCIIParamsTag = params;
end

%--------------------------------------------------------------------------

function pixelScale = constructModelPixelScale(R)
% Compute the GeoTIFF pixel scale, given a referencing object, R, with
% columns starting from the north. The standard GeoTIFF configuration is
% for columns to run from north-to-south (for a positive scaleY). This
% condition must be validated prior to calling this function.

if strcmp(R.CoordinateSystemType, 'geographic')
    % geographic
    if strcmp(R.RasterInterpretation, 'cells')
        % cells
        scaleX = R.CellExtentInLongitude;
        scaleY = R.CellExtentInLatitude;
    else
        % postings
        scaleX = R.SampleSpacingInLongitude;
        scaleY = R.SampleSpacingInLatitude;
    end
else
    % planar
    if strcmp(R.RasterInterpretation, 'cells')
        % cells
        scaleX = R.CellExtentInWorldX;
        scaleY = R.CellExtentInWorldY;
    else
        % postings
        scaleX = R.SampleSpacingInWorldX;
        scaleY = R.SampleSpacingInWorldY;
    end
end
pixelScale = [scaleX scaleY 0];

%--------------------------------------------------------------------------

function tiepoint = constructModelTiepoint(R)
% Compute the default GeoTIFF tie point.

% Map from the intrinsic system used in Mapping Toolbox to the GeoTIFF
% intrinsic system.
xiTie = 0;
yiTie = 0;
if strcmp(R.CoordinateSystemType, 'geographic')
    xwTie = R.LongitudeLimits(1);
    ywTie = R.LatitudeLimits(2);
else
    xwTie = R.XWorldLimits(1);
    ywTie = R.YWorldLimits(2);
end
tiepoint = [xiTie yiTie 0 xwTie ywTie 0];

%--------------------------------------------------------------------------

function [keys, S] = constructGeoKeyDirectory(GeoKeyTag)
% Construct the GeoKeyDirectory value and return it as an [M-by-4] double
% array in KEYS. S is a structure containing information to construct the
% GeoAsciiParamsTag and GeoDoubleParamsTag if needed.

% Sort the GeoKeyTag structure so that the GeoKey IDs are in ascending
% order.
[GeoKeyTag, GeoKeyID] = sortGeoKeyDirectoryTag(GeoKeyTag);

% Assign values for the key version information.
% Values are based on libgeotiff version 1.2.5
% All EPSG data is obtained from this library version.
% The values are found in the following headers:
%
%    Header file   Name               Value
%    ----------    ---                -----
%    geokeys.h    GvCurrentRevision     1
%    geovalues.h  GvCurrentMinorRev     0
%    geotiff.h    GvCurrentVersion      1
%    
KeyCurrentVersion = 1;
KeyRevisionMajor  = 1;
KeyRevisionMinor  = 0;

% Set numberOfKeys to number of GeoKeys in the GeoKeyID structure.
geoKeyNames = fieldnames(GeoKeyID);
numberOfKeys = numel(geoKeyNames);

% Create the M-by-4 keys array to hold all the GeoKeys and the header.
keys = zeros(size(numberOfKeys + 1, 4));

% Assign indices for access into KEY.
KeyID = 1;
KeyLocation = 2;
KeyCount = 3;
KeyValue = 4;

% Assign the header values.
keyNum = 1;
keys(keyNum, KeyID) = KeyCurrentVersion;
keys(keyNum, KeyLocation) = KeyRevisionMajor;
keys(keyNum, KeyCount) = KeyRevisionMinor;
keys(keyNum, KeyValue) = numberOfKeys;
keyNum = keyNum + 1;

% Initialize values.
asciiOffset  = 0;
doubleOffset = 0;
doubleIndex = false(size(geoKeyNames));
asciiIndex  = false(size(geoKeyNames));
GeoAsciiParamsTagID  = Tiff.TagID.GeoASCIIParamsTag;
GeoDoubleParamsTagID = Tiff.TagID.GeoDoubleParamsTag;
geoDoubleKeyNames = getGeoDoubleKeyNames;

% Loop through all the GeoKeys in the GeoKeyTag structure and assign
% the GeoKeys from the structure into the appropriate tag location.
for k=1:numel(geoKeyNames)
    keyName = geoKeyNames{k};
    keyID = GeoKeyID.(keyName);
    value = GeoKeyTag.(keyName);
    if ischar(value)
        % Use GeoAsciiParamsTag to store value.
        keyLocation = GeoAsciiParamsTagID;
        keyCount = numel(value) + 1;  % include '|' terminator
        keyValue = asciiOffset;
        asciiIndex(k) = true;
        asciiOffset = asciiOffset + keyCount;
    elseif isDoubleParam(keyName, value, geoDoubleKeyNames)
        % Use GeoDoubleParamsTag to store value.
        keyLocation = GeoDoubleParamsTagID;
        keyCount = numel(value);
        keyValue = doubleOffset;
        doubleIndex(k) = true;
        doubleOffset = doubleOffset + keyCount;        
    else
        % Use GeoKeyDirectoryTag to store value.
        keyLocation = 0;
        keyCount = 1;
        keyValue = value;
    end
    
    % Assign the GeoKey into the keys array.
    keys(keyNum, KeyID) = keyID;
    keys(keyNum, KeyLocation) = keyLocation;
    keys(keyNum, KeyCount) = keyCount;
    keys(keyNum, KeyValue) = keyValue;
    keyNum = keyNum + 1;
end

% Create a structure to return values for the GeoAsciiParamsTag and
% GeoDoubleParamsTag. Assign the cell array field value after structure
% construction to ensure a scalar structure.
S = struct( ...
    'DoubleOffset', doubleOffset, ...
    'DoubleNames',  '', ...
    'AsciiOffset',  asciiOffset, ...
    'AsciiNames',   '');

S.DoubleNames = geoKeyNames(doubleIndex);
S.AsciiNames  = geoKeyNames(asciiIndex);

%--------------------------------------------------------------------------

function keyNames = getGeoDoubleKeyNames()
% Return a cell array of all GeoDoubleParams key names as defined by the
% GeoTIFF specification at: 
% http://geotiff.maptools.org/spec/geotiff2.7.html#2.7

persistent geoDoubleParamsKeyNames
if isempty(geoDoubleParamsKeyNames) 
    
    geoDoubleParamsKeyNames = { ...
        'GeogPrimeMeridianLongGeoKey', ...
        'GeogLinearUnitSizeGeoKey', ...
        'GeogAngularUnitSizeGeoKey', ...
        'GeogSemiMajorAxisGeoKey', ...
        'GeogSemiMinorAxisGeoKey', ...
        'GeogInvFlatteningGeoKey', ...
        'ProjLinearUnitSizeGeoKey', ...
        'ProjStdParallel1GeoKey', ...
        'ProjStdParallel2GeoKey', ...
        'ProjLinearUnitSizeGeoKey', ...
        'ProjStdParallel1GeoKey', ...
        'ProjStdParallel2GeoKey', ...
        'ProjNatOriginLongGeoKey', ...
        'ProjNatOriginLatGeoKey', ...
        'ProjFalseEastingGeoKey', ...
        'ProjFalseNorthingGeoKey', ...
        'ProjFalseOriginLongGeoKey', ...
        'ProjFalseOriginLatGeoKey', ...
        'ProjFalseOriginEastingGeoKey', ...
        'ProjFalseOriginNorthingGeoKey', ...
        'ProjCenterLongGeoKey', ...
        'ProjCenterLatGeoKey', ...
        'ProjCenterEastingGeoKey', ...
        'ProjFalseOriginNorthingGeoKey', ...
        'ProjScaleAtNatOriginGeoKey', ...
        'ProjScaleAtCenterGeoKey', ...
        'ProjAzimuthAngleGeoKey', ...
        'ProjStraightVertPoleLongGeoKey '};
    geoDoubleParamsKeyNames = sort(geoDoubleParamsKeyNames);
end

keyNames = geoDoubleParamsKeyNames;
    
%--------------------------------------------------------------------------

function [GeoKeyTag, GeoKeyID] = sortGeoKeyDirectoryTag(GeoKeyTag)
% Sort the 'GeoKeyDirectoryTag' structure, GeoKeyTag, such that the field
% names of GeoKeyTag are in ascending order according to their ID number
% (to conform to the GeoTIFF specification for the GeoKeyDirectoryTag
% value). The GeoKeyID structure contains field names that match the field
% names in GeoKeyID. The values in GeoKeyID are the ID numbers. The result
% is such that the GeoKey ID of FIELDN of GeoKeyTag is the value of
% GeoKeyID.(FIELDN).

% Construct the GeoKey ID to name map.
idToNameMap = constructGeoKeyDirectoryMap;

% GeoKey names.
geoKeyNames = idToNameMap.values;

% Field names of GeoKeyTag
tagNames = fieldnames(GeoKeyTag);

% Index of matching GeoKey names to field names.
geoKeyIndex = ismember(geoKeyNames, tagNames);

% Sorted GeoKey ID numbers in the map. 
% Each idToNameMap key element corresponds to the idToNameMap values
% element. The map is created with ascending keys (GeoKey IDs); therefore
% the names are in GeoKeyID ascending order too.
geoKeyIDs = idToNameMap.keys;

% Sorted field names of GeoKeyTag.
sortedTagNames = geoKeyNames(geoKeyIndex);

% Sorted GeoKey ID numbers to match only the fields in GeoKeyTag.
sortedGeoKeyIDs = geoKeyIDs(geoKeyIndex);

% Sorted GeoKeyTag structure.
GeoKeyTag = orderfields(GeoKeyTag, sortedTagNames);

% Structure containing the GeoKey ID numbers as values for each field name
% in GeoKeyTag.
GeoKeyID = cell2struct(sortedGeoKeyIDs, sortedTagNames, 2);

%--------------------------------------------------------------------------

function tf = isDoubleParam(keyName, value, keyNames)
% Return true if VALUE requires storage in GeoDoubleParamsTag.
% VALUE requires double storage if it is:
%  1) not integer-valued, or
%  2) greater than the maximum value of an unsigned short int (uint16), or
%  3) less than the minimum value of an unsigned short int, or
%  4) not scalar
%  5) is a TYPE = DOUBLE geokey (defined from the GeoTIFF specification at
%       http://geotiff.maptools.org/spec/geotiff2.7.html#2.7)

% True if the value can fit within a C "short int".
isUint16 = isequal(double(value), double(uint16(value)));

% Must be a member of the cell array, keyNames (GeoDoubleParamsTag key),
% or not scalar-valued, or not fit within the size of a uint16 value. 
tf = any(strcmp(keyName, keyNames)) || ~isUint16 || ~isscalar(value);

%--------------------------------------------------------------------------

function params = constructGeoDoubleParams(GeoKey, geoKeyNames, numValues)
% Construct the corresponding parameters for the GeoKeyDoubleParamsTag by
% appending each value in GeoKey.(geoKeyNames{k}) to the output double 
% array, params.

params = zeros(1, numValues);
getValue = @(x)(GeoKey.(x));
params = appendKeyValues(params, getValue, geoKeyNames);

%--------------------------------------------------------------------------

function params = constructGeoAsciiParams(GeoKey, geoKeyNames, numValues)
% Construct the corresponding parameters for the GeoKeyAsciiParamsTag by
% appending each value in GeoKey.(geoKeyNames{k}) to the output, params.
% Each element value must be appended with the GeoTIFF termination
% character, '|'.

params(1:numValues) = ' ';
getValue = @(x)([GeoKey.(x) '|']);
params = appendKeyValues(params, getValue, geoKeyNames);

%--------------------------------------------------------------------------

function keyValues = appendKeyValues(keyValues, getValue, geoKeyNames)
% Append the values obtained from the function handle, getValue, to the
% keyValues parameter for each element in the cell array, geoKeyNames.

endIndex = 0;
for k=1:numel(geoKeyNames)
    keyName = geoKeyNames{k};
    value = getValue(keyName);
    startIndex = endIndex + 1;
    endIndex = startIndex + numel(value) - 1;
    keyValues(1,startIndex:endIndex) = value;
end

%--------------------------------------------------------------------------

function TiffTags = constructTiffTags(A, cmap, UserTags)
% Construct TIFF tags from R and the GeoKeyTag structure.

% Create TiffTags structure based on characteristics of input image.
TiffTags = assignTiffTags(A, cmap);

% Modify TiffTags based on user inputs.
TiffTags = updateTagsFromInputs(TiffTags, UserTags);

%--------------------------------------------------------------------------

function TiffTags = assignTiffTags(A, cmap)
% Construct the TiffTags structure based on image characteristics.

% Create structure of TIFF tags.
TiffTags = struct( ...
    'Photometric', [], ...
    'Compression', [], ...
    'ImageLength', [], ...
    'ImageWidth',  [], ...
    'SamplesPerPixel', [], ...
    'BitsPerSample',   [], ...
    'SampleFormat', [], ...
    'PlanarConfiguration', [], ...
    'Software', []); 

% Assign the tag value for Photometric.
if ~isempty(cmap)
    photometric = Tiff.Photometric.Palette;
    TiffTags.ColorMap = cmap;
else
    if ndims(A) == 3 && size(A,3) == 3 ...
            && any(strcmp(class(A), {'uint8', 'uint16'}))
        photometric = Tiff.Photometric.RGB;      
    else
        photometric = Tiff.Photometric.MinIsBlack;
    end
end
TiffTags.Photometric = photometric;            

% Assign the default tag value for Compression.
TiffTags.Compression = Tiff.Compression.PackBits;

% Assign the tag values for the image size.
TiffTags.ImageLength = size(A, 1);
TiffTags.ImageWidth  = size(A, 2);

% Assign the tag value for SamplesPerPixel.
map.internal.assert(ndims(A) <= 3, 'map:geotiff:invalidNumberOfDimensions');
TiffTags.SamplesPerPixel = size(A,3);

% Assign the tag values for BitsPerSample and SampleFormat.
classType = class(A);
TiffTags.BitsPerSample = assignBitsPerSample(classType);
TiffTags.SampleFormat  = assignSampleFormat(classType);

% Assign the tag value for PlanarConfiguration.
TiffTags.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

% Assign file organization tags.
TiffTags = assignFileOrganizationTags(TiffTags, A);

% Assign the tag value for Software.
v1 = ver('MATLAB');
v2 = ver('map');
software = [v1.Name ' ' v1.Version ', ' v2.Name, ' ' v2.Version];
TiffTags.Software = software;

%--------------------------------------------------------------------------

function bitsPerSample = assignBitsPerSample(classType)
% Assign the value for BitsPerSample based on class type.

switch classType
    case {'logical'}       
        bitsPerSample = 1;
        
    case {'uint8', 'int8'}
        bitsPerSample = 8;
        
    case {'uint16', 'int16'}
        bitsPerSample = 16;
        
    case {'uint32', 'int32', 'single'}
        bitsPerSample = 32;

    case {'double'}
        bitsPerSample = 64;

    otherwise
        error(message('map:geotiff:unknownClass', classType));
end

%--------------------------------------------------------------------------

function sampleFormat = assignSampleFormat(classType)
% Assign the value for SampleFormat based on class type.

switch classType
    case {'logical', 'uint8', 'uint16', 'uint32'}       
        sampleFormat = Tiff.SampleFormat.UInt;
                
    case {'int8', 'int16', 'int32'}
        sampleFormat = Tiff.SampleFormat.Int;
                    
    case {'single', 'double'}
        sampleFormat = Tiff.SampleFormat.IEEEFP;

    otherwise
        error(message('map:geotiff:unknownClass', classType));
end

%--------------------------------------------------------------------------

function TiffTags = assignFileOrganizationTags(TiffTags, A)
% Assign the file organization tags, TileLength and TileWidth or
% RowsPerStrip.

% Set the tile size if the image is at least 160-by-160. Note the minimum
% tile dimension is 16-by-16. TileLength and TileWidth must be multiples of
% 16. Create tiles such that the maximum number of tiles is 100 (which can
% be obtained from the Tiff class method, numberOfTiles).
if size(A,1) >= 160 && size(A,2) >= 160
    maxTileWidth  = size(A,2)/10;
    maxTileLength = size(A,1)/10;
    tileWidth  = fix(maxTileWidth/16)*16;
    tileLength = fix(maxTileLength/16)*16;
    
    if tileWidth ~= maxTileWidth
        % Increase tileWidth to next multiple of 16.
        tileWidth = tileWidth + 16;
    end
    if tileLength ~= maxTileLength
        % Increase tileLength to next multiple of 16.
        tileLength = tileLength + 16;
    end
   
    TiffTags.TileWidth  = tileWidth;
    TiffTags.TileLength = tileLength;
else
    TiffTags.RowsPerStrip = 1;
end

%--------------------------------------------------------------------------

function TiffTags = updateTagsFromInputs(TiffTags, UserTags)
% Copy the fields from UserTags to TiffTags.

userNames = fieldnames(UserTags);

% Remove TileLength and TileWidth if RowsPerStrip has been specified.
tileNames = {'TileLength', 'TileWidth'};
if any(strcmpi('RowsPerStrip', userNames)) ...
        && all(isfield(TiffTags, tileNames))
   TiffTags = rmfield(TiffTags, tileNames);
end

% Remove Software if the user requested the empty string.
field = 'Software';
index = strcmpi(field, userNames);
if any(index) && isempty(UserTags.(field))
    TiffTags = rmfield(TiffTags, field);
    UserTags = rmfield(UserTags, field);
    userNames(index) = [];
end

% Copy or create the fields from UserTags to TiffTags.
for k=1:numel(userNames)
    field = userNames{k};
    TiffTags.(field) = UserTags.(field);
end


%--------------------------------------------------------------------------

function writeGeoTiffFile(filename, A, TiffTags, GeoTiffTags, RPCCoefficientTag)

% Open output TIFF file.
t = Tiff(filename, 'w');

% In case of error, or when terminating the function, close the TIFF file.
tobj = onCleanup(@()close(t));

% Set the tags and write image data.
t.setTag(TiffTags);
t.setTag(GeoTiffTags);
t.write(A);

% Close the Tiff object.
delete(tobj)

% Add RPCCoefficientTag if not empty.
if ~isempty(RPCCoefficientTag)
    tagValue = double(RPCCoefficientTag);
    map.geotiff.internal.writeRPCCoefficientTag(filename, tagValue);
end
