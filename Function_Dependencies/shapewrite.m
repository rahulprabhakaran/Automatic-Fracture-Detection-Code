function shapewrite(varargin)
%SHAPEWRITE Write geographic vector data to shapefile
%
%   SHAPEWRITE(S, FILENAME) writes the vector geographic features stored in
%   S to disk in shapefile format.  S is either a mappoint vector, mapshape
%   vector, mapstruct (with 'X' and 'Y' coordinate fields), geopoint
%   vector, geoshape vector, or a geostruct (with 'Lat' and 'Lon' fields)
%   with the following restrictions on its attribute fields:
%
%   * Each attribute field value must be either a real, finite, scalar
%     double or a character string.  
%
%   * The type of a given attribute must be consistent across all features.
%
%   If S is a geopoint vector, geoshape vector, or a geostruct, then the
%   latitude and longitude values are written out as 'Y' and 'X'
%   coordinates, respectively, in the shapefile.
%
%   FILENAME must be a character string specifying the output
%   file name and location.  If an extension is included, it must be
%   '.shp' or '.SHP'. SHAPEWRITE creates three output files,
%
%                        [BASENAME '.shp']
%                        [BASENAME '.shx']
%                        [BASENAME '.dbf']
%
%   where BASENAME is FILENAME without its extension.
%
%   If a given attribute is integer-valued for all features, then it is
%   written to the [BASENAME '.dbf'] file as an integer.  If an attribute
%   is a non-integer for any feature, then it is written as a fixed point
%   decimal value with six digits to the right of the decimal place. 
%
%   SHAPEWRITE(S, FILENAME, 'DbfSpec', DBFSPEC) writes a shapefile in which
%   the content and layout of the DBF file is controlled by a DBF
%   specification, indicated here by the parameter value DBFSPEC.  A DBF
%   specification is a scalar MATLAB structure with one field for each
%   feature attribute to be included in the output shapefile.  To include
%   an attribute in the output, make sure to provide a field in DBFSPEC
%   with a fieldname identical to the attribute name (the corresponding
%   fieldname in S), and assign to that field a scalar structure with the
%   following four fields:
% 
%     FieldName -- The field name to be used in the file
%
%     FieldType -- The field type to be used in the file ('N' or 'C')
%
%     FieldLength -- The field length in the file, in bytes
%
%     FieldDecimalCount -- For numeric fields, the number of digits to the
%                          right of the decimal place
%
%   When a DBF spec is provided, a given attribute will be included in the
%   output file only if it matches the name of a field in the spec.
%
%   The easiest way to construct a DBF spec is to call MAKEDBFSPEC, then
%   modify the output to remove attributes or change the FieldName,
%   FieldLength, or FieldDecimalCount for one or more attributes.  See the
%   help for MAKEDBFSPEC for more details and an example.  
%
%   Example
%   -------
%   % Derive a shapefile from concord_roads.shp in which roads of CLASS 5
%   % and greater are omitted.  Note the use of the 'Selector' option in
%   % shaperead, together with an anonymous function, to read only the main
%   % roads from the original shapefile.
%   shapeinfo('concord_roads')  % 609 features
%   S = shaperead('concord_roads', 'Selector', ...
%                 {@(roadclass) roadclass < 4, 'CLASS'});
%   shapewrite(S, 'main_concord_roads.shp')
%   shapeinfo('main_concord_roads')  % 107 features
%
%   See also MAKEDBFSPEC, SHAPEINFO, SHAPEREAD, UPDATEGEOSTRUCT.

% Copyright 2003-2016 The MathWorks, Inc.

narginchk(2, Inf);

[S, basename, dbfspec] = parseInputs(varargin{:});
[shapeType, boundingBox, index] = writeSHP(S,basename);
writeSHX(shapeType, boundingBox, index, basename);
if ~isempty(dbfspec)
    dbfwrite(S, basename, dbfspec)
end
 
%--------------------------------------------------------------------------

function [shapeType, boundingBox, index] = writeSHP(S,basename)
% Write the main (SHP) file.

% Open the SHP file.
fid = fopen([basename '.shp'],'w','ieee-be');
if fid < 0
    error(message('map:shapefile:failedToOpenFile', [ basename, '.shp' ]))
end

% Get the shape type and a handle to function to write individual records.
[shapeType, writefcn] = getShapeType(S);

% Determine coordinate fields.
[xField, yField] = coordinateFieldNames(S);

% Write 100 bytes of zeros to reserve room for the file header.
fwrite(fid, uint8(zeros(1,100)), 'uint8');

% Write an SHP record for each element in S.
% Accumulate an index array and bounding box.
headerLengthInWords = 50;
fileLengthInWords = headerLengthInWords;
boundingBox = [Inf Inf -Inf -Inf];
index = zeros(2,length(S));
for k = 1:length(S)
    x = S(k).(xField);
    y = S(k).(yField);
    [fileLengthInWords, boundingBox, index(:,k)] = ...
        shpWriteRecord(fid, writefcn, shapeType, x, y, k, fileLengthInWords, boundingBox);
end

% Back up to the beginning of the file and write the header into the first 100 bytes.
fseek(fid, 0, 'bof');
shpWriteHeader(fid,shapeType,fileLengthInWords,boundingBox);

% Close the SHP file.
fclose(fid);

%--------------------------------------------------------------------------

function [xField, yField] = coordinateFieldNames(S)
% Determine coordinate field / property names from input struct / class

if isstruct(S)
    if isfield(S,'X') && isfield(S,'Y')
        % mapstruct input
        xField = 'X';
        yField = 'Y';
    elseif isfield(S,'Lon') && isfield(S,'Lat')
        % geostruct input
        xField = 'Lon';
        yField = 'Lat';
    else
        error(message('map:geostruct:missingCoordinateFields'))
    end
else
    classname = class(S);
    if any(strcmp(classname,{'mappoint','mapshape'}))
        % mappoint or mapshape input
        xField = 'X';
        yField = 'Y';
    else
        % geopoint or geoshape input
        xField = 'Longitude';
        yField = 'Latitude';
    end
end

%--------------------------------------------------------------------------

function [fileLengthInWords, boundingBox, index] = shpWriteRecord(fid,...
      writefcn, shapeType, x, y, recordNumber, fileLengthInWords, boundingBox)
% Write an SHP record for each structure element s.  Return index entry.

recordHeaderLengthInWords = 4;
offsetInWords = ftell(fid) / 2;
[contentLengthInWords, recordBoundingBox] ...
    = writefcn(fid, recordNumber, shapeType, x, y);
fileLengthInWords = fileLengthInWords ...
    + recordHeaderLengthInWords + contentLengthInWords;
boundingBox(1:2) = min([boundingBox(1:2); recordBoundingBox(1:2)]);
boundingBox(3:4) = max([boundingBox(3:4); recordBoundingBox(3:4)]);
index = [offsetInWords; contentLengthInWords];

%--------------------------------------------------------------------------

function writeSHX(shapeType, boundingBox, index, basename)
% Write the SHX file.

% Open the SHX file.
filename = [basename '.shx'];
fid = fopen(filename,'w','ieee-be');
if fid < 0
    error(message('map:shapefile:failedToOpenFile', filename))
end

% Write the header.
headerLengthInWords = 50;
fileLengthInWords = headerLengthInWords + 4 * size(index,2);
shpWriteHeader(fid,shapeType,fileLengthInWords,boundingBox);

% Write the index records.
fwrite(fid,int32(index),'int32','ieee-be');

% Close the SHX file.
fclose(fid);

%--------------------------------------------------------------------------

function [shapeType, writefcn] = getShapeType(S)
% Return the numerical shape type code to be written to the .SHP file and
% the handle to a function for writing an individual .SHP record.

if isstruct(S)
   % S is mapstruct or geostruct
    switch(lower(S(1).Geometry))
        case 'point'
            shapeType = 1;
            writefcn = @shpWritePoint;
            
        case 'multipoint'
            shapeType = 8;
            writefcn = @shpWriteMultiPoint;
            
        case {'line','polyline'}
            shapeType = 3;
            writefcn = @shpWritePoly;
            
        case 'polygon'
            shapeType = 5;
            writefcn = @shpWritePoly;
    end
elseif any(strcmp(class(S), {'mappoint','geopoint'}))
    shapeType = 1;
    writefcn = @shpWritePoint;
else
    % S is a mapshape or geoshape
    switch(S.Geometry)
        case 'point'
            shapeType = 8;
            writefcn = @shpWriteMultiPoint;
            
        case 'line'
            shapeType = 3;
            writefcn = @shpWritePoly;
            
        case 'polygon'
            shapeType = 5;
            writefcn = @shpWritePoly;
    end
end

%--------------------------------------------------------------------------

function shpWriteHeader(fid, shapeType, fileLengthInWords, boundingBox)

fileCode = 9994;
version  = 1000;

bytes0thru27 = int32([fileCode 0 0 0 0 0 fileLengthInWords]);
fwrite(fid, bytes0thru27, 'int32', 'ieee-be');
% Check that count = 7

bytes28thru35 = int32([version shapeType]);
fwrite(fid, bytes28thru35, 'int32', 'ieee-le');
% Check that count = 2

bytes36thru99 = [boundingBox 0 0 0 0];
fwrite(fid, bytes36thru99, 'double', 'ieee-le');
% Check that count = 8

%--------------------------------------------------------------------------

function [contentLengthInWords, boundingBox] ...
    = shpWritePoint(fid,recordNumber,shapeType,X,Y)
% Write a Point record. Return the record content length measured in
% 16-bit words and the bounding box for the current record.

% Assign a degenerate bounding box for use in calculating the bounding box
% for the overall file.
boundingBox = [X Y X Y]; 

% The content length (the length of the record excluding its header) is
% fixed at 20 bytes (10 words) for Point shapes.
contentLengthInWords = 10;

% Write the record header
fwrite(fid, int32([recordNumber contentLengthInWords]), 'int32', 'ieee-be');

% Write the shape type.
fwrite(fid, int32(shapeType), 'int32', 'ieee-le');

% Write the point coordinates.
fwrite(fid, [X Y], 'double', 'ieee-le');

%--------------------------------------------------------------------------

function [contentLengthInWords, boundingBox] ...
    = shpWriteMultiPoint(fid,recordNumber,shapeType,X,Y)
% Write a MultiPoint record. Return the record content length measured in
% 16-bit words and the bounding box for the current record.

% Calculate the 2-by-numPoints points array and the bounding box.
points = [X(:)'; Y(:)'];
boundingBox = [min(X) min(Y) max(X) max(Y)]; 

% Calculate the content length (the length of the record excluding its header).
contentLengthInWords = (40 + 8 * numel(points)) / 2;

% Write the record header
fwrite(fid, int32([recordNumber contentLengthInWords]), 'int32', 'ieee-be');

% Write the shape type.
fwrite(fid, int32(shapeType), 'int32', 'ieee-le');

% Write the bounding box.
fwrite(fid, boundingBox, 'double', 'ieee-le');

% Write the number of points.
fwrite(fid, size(points,2), 'int32', 'ieee-le');

% Write the point coordinates.
fwrite(fid, points, 'double', 'ieee-le');

%--------------------------------------------------------------------------

function [contentLengthInWords, boundingBox] ...
    = shpWritePoly(fid,recordNumber,shapeType,X,Y)
% Write a PolyLine or Polygon record. Return the record content length measured in
% 16-bit words and the bounding box for the current record.

% Calculate the 1-by-numParts parts index, the 2-by-numPoints points array,
% and the bounding box.
[parts, points] = constructPartsIndexAndPointsArray(X,Y);
boundingBox = [min(X) min(Y) max(X) max(Y)]; 

% Calculate the content length (the length of the record excluding its header).
contentLengthInWords = (44 + 4 * numel(parts) + 8 * numel(points)) / 2;

% Write the record header.
fwrite(fid, int32([recordNumber contentLengthInWords]), 'int32', 'ieee-be');

% Write the shape type.
fwrite(fid, int32(shapeType), 'int32', 'ieee-le');

% Write the bounding box.
fwrite(fid, boundingBox, 'double', 'ieee-le');

% Write the part and point counts, and the parts index.
fwrite(fid, int32([numel(parts) size(points,2) parts]), 'int32','ieee-le');

% Write the point coordinates.
fwrite(fid, points, 'double', 'ieee-le');

%--------------------------------------------------------------------------

function [parts, points] = constructPartsIndexAndPointsArray(X,Y)
% Calculate the zero-based index array PARTS giving the offset to the where
% the start of each non-NaN sequence in vector X will be after all the NaNs
% were removed from X.  Organize the non-NaN elements of X and Y into a
% 2-by-numPoints array POINTS with X-coordinates in row 1 and Y-coordinates
% in row 2.

% Check that isequal(isnan(X),isnan(Y)). If not, then throw an error here
% and catch it above so we have a chance to close open files.  Also require
% that X and Y each have at least one non-NaN element.

% Determine where the NaNs are and add a (virtual) terminating NaN.
nanplaces = [find(isnan(X(:))); 1+numel(X)];

% Compute the parts index. Note that in both the initialization of parts and
% in subsequent updates, we calculate the offset to the next part before
% determining if that part exists.  This is fine, except that at the end we
% end up with parts(end) containing the offset of a non-existent part, so
% we must remove it.
parts = 0;
if any(nanplaces < numel(X))
    indexOfLastNaN = 0;
    for k = 1:numel(nanplaces)
        partLength = nanplaces(k) - (indexOfLastNaN + 1);
        if partLength > 0
            parts(end+1) = parts(end) + partLength; %#ok<AGROW>
        end
        % If partLength is zero, then the current NaN was preceded
        % immediately by another NaN, and should hence be ignored.
        indexOfLastNaN = nanplaces(k);
    end
    parts(end) = [];
end

% Calculate the points array:
%    [X(1) X(2) ...;
%     Y(1) Y(2) ...]
% after removing NaNs.
x = X(:)'; x(isnan(x)) = [];
y = Y(:)'; y(isnan(y)) = [];
points = [x; y];

%--------------------------------------------------------------------------

function dbfwrite(S, basename, dbfspec)
% Write the DBF file.

% Open the DBF file.
filename = [basename '.dbf'];
fid = fopen(filename,'w','ieee-le');
if fid < 0
    error(message('map:shapefile:failedToOpenFile', filename))
end

% Determine the record length, construct format string, and for each field
% in the dbfspec, determine its position among the fields of the geographic
% data structure array, S.
[S, recordLength, fmt, geostructFieldIndex, nanmask] ...
    = setupRecordLayout(dbfspec, S);

% Write the DBF file: header, data records, and EOF marker
numberOfFeatures = length(S(:));
dbfWriteTableFileHeader(fid, numberOfFeatures, dbfspec, recordLength);
dbfWriteAllRecords(fid, S, geostructFieldIndex, fmt, nanmask)
dbfMarkEOF(fid);

% Close the DBF file.
fclose(fid);

%--------------------------------------------------------------------------

function dbfWriteAllRecords(fid, S, geostructFieldIndex, fmt, nanmask)

numberOfFeatures = length(S(:));
if isstruct(S)
    for k = 1:numberOfFeatures
        allFields = struct2cell(S(k));
        values = allFields(geostructFieldIndex);
        dbfWriteTableRecord(fid, values, fmt, nanmask);
    end
else
    % Convert the dynamic vector to a scalar structure.
    S = struct(S);
    values = cell(length(geostructFieldIndex), 1);
    p = fieldnames(S);
    p = p(geostructFieldIndex);
    
    % Loop through each feature and convert the requested attributes
    % (defined by geostructFieldIndex) to values in a cell array. Each
    % attribute is either a scalar numeric value, a value in a cell array,
    % or a string value.
    for k = 1:numberOfFeatures
        for n = 1:length(geostructFieldIndex)
            feature = S.(p{n});
            if iscell(feature)
                values{n, 1} = feature{k};
            elseif ischar(feature)
                values{n, 1} = feature;
            else
                values{n, 1} = feature(k);
            end
        end
        dbfWriteTableRecord(fid, values, fmt, nanmask);
    end
end

%--------------------------------------------------------------------------

function [S, recordLength, fmt, geostructFieldIndex, nanmask] = ...
    setupRecordLayout(dbfspec, S)

% While iterating to find the field length and building up a format
% specification string one attribute at a time, also construct a cell array,
% nulldata, which contains the value NaN for each numeric-valued
% attribute and the value '' for each string-valued attribute.
% This function modifies the string-valued attributes in the input
% geostruct array, S. For each character type attribute, convert the
% Unicode characters to their native representation, using the default
% character encoding stream.

if isstruct(S)
    geostructFields = fields(S);
else
    geostructFields = properties(S);
end

recordLength = 1;  % Account for leading space
dbfFieldSpecs = struct2cell(dbfspec);
fmt = cell(size(dbfFieldSpecs));
nulldata = fmt;
attributeNames = fields(dbfspec);
geostructFieldIndex = zeros(size(attributeNames));
for k = 1:numel(dbfFieldSpecs)
    f = dbfFieldSpecs{k};
    switch(f.FieldType)
        case 'N'
            fmt{k} = sprintf('%s%d.%df', '%', f.FieldLength, f.FieldDecimalCount);
            nulldata{k} = NaN;
        otherwise % 'C'
            fmt{k} = sprintf('%s%ds', '%-', f.FieldLength);
            nulldata{k} = '';
            S = convertUnicodeToDefault(S, attributeNames{k});
    end
    recordLength = recordLength + f.FieldLength;
    geostructFieldIndex(k) = ...
        find(strcmp(attributeNames{k}, geostructFields));
end

% Concatenate the format strings for the individual fields into a single
% format string for the entire record, including a leading space.
fmt = [' ', fmt{:}];

% nanmask is a string, with length equal to recordLength.  It is used
% in subfunction dbfWriteTableRecord to ensure that each instance of a
% NaN-valued attribute is represented by blanks (space characters) in the
% xBase (.DBF) file.  nanmask consists of substrings equal to 'NaN'
% interspersed among a background of space characters (' ').  There is
% one such substring for each numerical attribute.  Each 'NaN' substring
% located within nanmask at exactly the location where 'NaN' would be
% written by a call to sprintf(fmt,c) if c were a cell array with a
% numerical value of NaN for the corresponding attribute field.  nanmask
% might look something like this:  '     NaN  NaN        NaN    ', for
% example, assuming three numeric-valued attributes.
nanmask = sprintf(fmt, nulldata{:});

%--------------------------------------------------------------------------

function S = convertUnicodeToDefault(S, attributeName)
% Convert the Unicode characters to the default encoding scheme.

if isstruct(S)
    for k = 1:numel(S)
        S(k).(attributeName) = ...
            char(unicode2native(S(k).(attributeName)));
    end
else
    values = S.(attributeName);
    if iscell(values)
        for k = 1:length(values)
            values{k} = ...
                char(unicode2native(values{k}));
        end
        S.(attributeName) = values;
    else
        S.(attributeName) = char(unicode2native(values));
    end
end

%--------------------------------------------------------------------------

function dbfWriteTableFileHeader(fid, numberOfRecords, dbfspec, recordLength)

version = 3;

dateVector = datevec(date);
year  = dateVector(1);
month = dateVector(2);
day   = dateVector(3);

dbfFieldSpecs = struct2cell(dbfspec);
nFields = numel(dbfFieldSpecs);
headerLength = 32 + 32 * nFields + 1;

% Bytes 0-3: dBASE version and today's date
fwrite(fid, [version year-1900 month day], 'uint8');

% Bytes 4-7: Number of records in the table
fwrite(fid, numberOfRecords, 'uint32');

% Bytes 8-9: Number of bytes in the header
fwrite(fid, headerLength, 'uint16');

% Bytes 10-11: Number of bytes per record
fwrite(fid, recordLength, 'uint16');

% Bytes 12-31: Reserved (zero-fill)
fwrite(fid, zeros(1,20), 'uint8');

% Bytes 32-n: Table field descriptors (n = 32 + 32 * nFields)
for k = 1:nFields
    dbfWriteTableFieldDescriptor(fid, dbfFieldSpecs{k});
end

% Byte n+1: Field terminator
fwrite(fid, hex2dec('0D') ,'uint8');

%--------------------------------------------------------------------------

function dbfWriteTableFieldDescriptor(fid, dbfFieldSpec)

% Bytes 0-10: Field name in ASCII
%   (truncate or zero-fill to precisely 11 bytes)
fwrite(fid, truncateOrFill(dbfFieldSpec.FieldName,11), 'uchar');

% Byte 11: Field type in ASCII ('C' or 'N')
fwrite(fid, dbfFieldSpec.FieldType, 'uchar');

% Bytes 12-15: Field data address (zero-fill)
fwrite(fid, zeros(1,4), 'uint8');

% Byte 16: Field length
fwrite(fid, dbfFieldSpec.FieldLength, 'uint8');

% Byte 17: Field decimal count
%   (number of digits to the right of the decimal place)
fwrite(fid, dbfFieldSpec.FieldDecimalCount, 'uint8');

% Bytes 18-31: Reserved or zero
%   ("work area ID, SET FIELDS flag", ".MDX field flag")
fwrite(fid, zeros(1,14), 'uint8');

%--------------------------------------------------------------------------

function str = truncateOrFill(str,len)
% Ensure a string containing precisely LEN characters,
% padded with zeros if necessary.

str(len+1:end) = [];
str(end+1:len) = 0;

%--------------------------------------------------------------------------

function dbfWriteTableRecord(fid, s, fmt, nanmask)
% Write a single record, including the preceding SPACE character.
% s is a scalar geostruct or mapstruct, geostructFieldIndex is an array of
% indices indicating which fields of s need to be written, and
% fmt is a format string with one element per field in
% geostructFieldIndex.  See the comment in subfunction setupRecordLayout
% for a definition of nanmask.

% Note that if s contains any NaN-valued numerical attributes, then the
% call to sprintf will result in the substring 'NaN' being written into
% the record string.  Any such substrings must be replaced with blanks
% (space characters) before writing the record to the xBase (.DBF) file.
% This is accomplished via logical indexing with the nanmask string.

record = sprintf(fmt, s{:,1});
record(record == nanmask) = ' ';
fwrite(fid, record, 'uchar');

%--------------------------------------------------------------------------

function dbfMarkEOF(fid)
% Mark end of dBASE file.

eofMarker = hex2dec('1A');
fwrite(fid, eofMarker, 'uint8');

%--------------------------------------------------------------------------

function [S, basename, dbfspec] = parseInputs(varargin)

% S is the first of two required inputs.
S = varargin{1};
validateattributes(S, ...
    {'struct', 'mappoint', 'geopoint', 'mapshape', 'geoshape'}, ...
    {'nonempty', 'vector'}, mfilename, 'S', 1)

% FILENAME is the second of two required inputs.
filename = varargin{2};
basename = validateFilename(filename);

% Identify and validate the parameter name-value pairs beginning with
% the third argument.
validParameterNames = {'DbfSpec'};
for k = 3:2:nargin
    parName = validatestring(varargin{k}, validParameterNames, ...
        mfilename, sprintf('PARAM%d',(k-1)/2), k); 
    checkExistence(k, nargin, 'a scalar structure', parName);
    dbfspec = varargin{k+1};
end

% If the user provided a DBF spec, validate it.  Otherwise, derive one. 
%
% Note: The makedbfspec function is called in both branches
%       (validateDbfSpec always invokes it), and it validates the geometry
%       and the attribute fields of S.
%
if exist('dbfspec','var')
    % Validate dbfspec and S and ensure consistency.
    dbfspec = validateDbfSpec(dbfspec, S);
else
    dbfspec = makedbfspec(S);
end

%--------------------------------------------------------------------------

function checkExistence(position, nargs, propertyDescription, propertyName)
% Error if missing the property value following a property name.

if (position + 1 > nargs)
    error(message('map:shapefile:missingParameterValue', propertyDescription, propertyName));
end

%--------------------------------------------------------------------------

function basename = validateFilename(filename)

validateattributes(filename, {'char'}, {'vector'}, mfilename, 'FILENAME', 2);
[pathstr, name, ext] = fileparts(filename);
if ~any(strcmpi(ext,{'','.shp'}))
    error(message('map:shapefile:invalidShpExtension', mfilename, filename))
end
basename = fullfile(pathstr,name);

%--------------------------------------------------------------------------

function dbfspec = validateDbfSpec(dbfspec, S)

% If dbfspec is empty, make sure it's an empty struct.
if isempty(dbfspec)
    dbfspec = struct([]);
    return  % No need to check anything else, including the attribute fields of S.
end

% Make sure that dbfspec is a structure.
map.internal.assert(isstruct(dbfspec), 'map:shapefile:dbfspecNotAStructure')

% Make sure that dbfspec is a scalar (or empty).
if numel(dbfspec) > 1
    error(message('map:shapefile:nonScalarDbfSpec'))
end

% Validate S, then make sure that dbfspec and S are mutually consistent.
defaultspec = makedbfspec(S);  % Validates attribute values in S
attributeNamesInS = fields(defaultspec);
attributeNamesInDbfSpec = fields(dbfspec);

% Check 1:  Every attribute in dbfspec must be present in S.
missingAttributes = setdiff(attributeNamesInDbfSpec,attributeNamesInS);
map.internal.assert(isempty(missingAttributes), 'map:shapefile:missingAttributes')

% Check 2:  Field types in dbfspec must be consistent with S.
%   While in this loop, use the default to fill in any fields missing from
%   dbfspec and ensure that the field length is at least 2 for
%   character-valued attributes and 3 for numerical attributes.
minCharFieldLength = 2;
minNumericalFieldLength = 3;  % Large enough to hold 'NaN'
for k = 1:numel(attributeNamesInDbfSpec)
    attributeName = attributeNamesInDbfSpec{k};
    fSpec = dbfspec.(attributeName);
    map.internal.assert(isstruct(fSpec) && isscalar(fSpec), ...
        'map:shapefile:expectedScalarDbfStructure', attributeName);    
    fDefault = defaultspec.(attributeName);
    
    if ~isfield(fSpec,'FieldName')
        dbfspec.(attributeName).FieldName = fDefault.FieldName;
    end

    if ~isfield(fSpec,'FieldType')
        dbfspec.(attributeName).FieldType = fDefault.FieldType;
        fSpec.FieldType = fDefault.FieldType;
    elseif fSpec.FieldType ~= fDefault.FieldType
        error(message('map:shapefile:extraneousAttributes', attributeName))
    end
    
    if ~isfield(fSpec,'FieldLength')
        dbfspec.(attributeName).FieldLength = fDefault.FieldLength;
    elseif fSpec.FieldType == 'C'
        dbfspec.(attributeName).FieldLength ...
            = max(minCharFieldLength, fSpec.FieldLength);
    else
        dbfspec.(attributeName).FieldLength ...
            = max(minNumericalFieldLength, fSpec.FieldLength);
    end
    
    if ~isfield(fSpec,'FieldDecimalCount')
        dbfspec.(attributeName).FieldDecimalCount = fDefault.FieldDecimalCount;
    end
end
