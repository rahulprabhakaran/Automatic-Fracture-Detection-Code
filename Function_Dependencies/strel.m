%STREL Create morphological structuring element.
%   SE = STREL('arbitrary',NHOOD) creates a flat structuring element with
%   the specified neigbhorhood.  NHOOD is a matrix containing 1's and 0's;
%   the location of the 1's defines the neighborhood for the morphological
%   operation.  The center (or origin) of NHOOD is its center element,
%   given by FLOOR((SIZE(NHOOD) + 1)/2).  You can also omit the 'arbitrary'
%   string and just use STREL(NHOOD).
%
%   SE = STREL('diamond',R) creates a flat diamond-shaped structuring
%   element with the specified size, R.  R is the distance from the
%   structuring element origin to the points of the diamond.  R must be a
%   nonnegative integer scalar.
%
%   SE = STREL('disk',R,N) creates a flat disk-shaped structuring element
%   with the specified radius, R.  R must be a nonnegative integer.  N must
%   be 0, 4, 6, or 8.  When N is greater than 0, the disk-shaped structuring
%   element is approximated by a sequence of N (or sometimes N+2)
%   periodic-line structuring elements.  When N is 0, no approximation is
%   used, and the structuring element members comprise all pixels whose
%   centers are no greater than R away from the origin.  N can be omitted,
%   in which case its default value is 4.  Note: Morphological operations
%   using disk approximations (N>0) run much faster than when N=0.  Also,
%   the structuring elements resulting from choosing N>0 are suitable for
%   computing granulometries, which is not the case for N=0.  Sometimes it
%   is necessary for STREL to use two extra line structuring elements in the
%   approximation, in which case the number of decomposed structuring
%   elements used is N+2.
%
%   SE = STREL('line',LEN,DEG) creates a flat linear structuring element
%   that is symmetric with respect to the neighborhood center. DEG
%   specifies the angle (in degrees) of the line as measured in a
%   counterclockwise direction from the horizontal axis. LEN is
%   approximately the distance between the centers of the structuring
%   element members at opposite ends of the line.
%
%   SE = STREL('octagon',R) creates a flat octagonal structuring element
%   with the specified size, R.  R is the distance from the structuring
%   element origin to the sides of the octagon, as measured along the
%   horizontal and vertical axes.  R must be a nonnegative multiple of 3.
%
%   SE = STREL('rectangle',MN) creates a flat rectangle-shaped structuring
%   element with the specified size.  MN must be a two-element vector of
%   nonnegative integers.  The first element of MN is the number rows in the
%   structuring element neighborhood; the second element is the number of
%   columns.
%
%   SE = STREL('square',W) creates a square structuring element whose
%   width is W pixels.  W must be a nonnegative integer scalar.
%
%   SE = STREL('sphere',R) creates a sphere-shaped structuring element
%   whose radius is R pixels. R must be a nonnegative integer scalar.
%
%   SE = STREL('cube',W) creates a cube structuring element whose
%   width is W pixels. W must be a nonnegative integer scalar.
%
%   SE = STREL('cuboid', XYZ) creates a cuboidal structuring element of
%   size XYZ. XYZ must be a three-element vector of nonnegative integers
%   specifying the number of rows, number of columns, and the number of
%   planes in the third dimension.
%
%   Notes
%   -----
%   For all shapes except 'arbitrary', structuring elements are constructed
%   using a family of techniques known collectively as "structuring element
%   decomposition."  The principle is that dilation by some large
%   structuring elements can be computed faster by dilation with a sequence
%   of smaller structuring elements.  For example, dilation by an 11-by-11
%   square structuring element can be accomplished by dilating first with a
%   1-by-11 structuring element and then with an 11-by-1 structuring
%   element.  This results in a theoretical performance improvement of a
%   factor of 5.5, although in practice the actual performance improvement
%   is somewhat less.  Structuring element decompositions used for the
%   'disk' and 'ball' shapes are approximations; all other decompositions
%   are exact.
%
%   The following syntaxes will still work, but are not recommended for
%   use. The non-flat structuring element shapes listed below can be
%   created using OFFSETSTREL.
%
%   SE = STREL('pair',OFFSET) creates a flat structuring element containing
%   two members.  One member is located at the origin; the second member's
%   location is specified by the vector OFFSET.  OFFSET must be a
%   two-element vector of integers.
%
%   SE = STREL('periodicline',P,V) creates a flat structuring element
%   containing 2*P+1 members.  V is a two-element vector containing
%   integer-valued row and column offsets.  One structuring element member
%   is located at the origin.  The other members are located at 1*V, -1*V,
%   2*V, -2*V, ..., P*V, -P*V.
%
%   SE = STREL('arbitrary',NHOOD,HEIGHT) creates a nonflat structuring
%   element with the specified neighborhood.  HEIGHT is a matrix the same
%   size as NHOOD containing the height values associated with each nonzero
%   element of NHOOD.  HEIGHT must be real and finite-valued.  You can also
%   omit the 'arbitrary' string and just use STREL(NHOOD,HEIGHT).
%
%   SE = STREL('ball',R,H,N) creates a nonflat "ball-shaped" (actually an
%   ellipsoid) structuring element whose radius in the X-Y plane is R and
%   whose height is H.  R must be a nonnegative integer, and H must be a
%   real scalar.  N must be an even nonnegative integer.  When N is greater
%   than 0, the ball-shaped structuring element is approximated by a
%   sequence of N nonflat line-shaped structuring elements.  When N is 0, no
%   approximation is used, and the structuring element members comprise all
%   pixels whose centers are no greater than R away from the origin, and the
%   corresponding height values are determined from the formula of the
%   ellipsoid specified by R and H.  If N is not specified, the default
%   value is 8.  Note: Morphological operations using ball approximations
%   (N>0) run much faster than when N=0.
%
%   Examples
%   --------
%       se1 = strel('square', 11)      % 11-by-11 square
%       se2 = strel('line', 10, 45)     % length 10, angle: 45 degrees
%
%       % Create and display a 2-D disk-shaped structuring element.
%       se3 = strel('disk', 15)
%       figure, imshow(se3.Neighborhood)
%
%       % Create and display a 3-D sphere-shaped structuring element.
%       se4 = strel('sphere', 15)
%       figure, isosurface(se4.Neighborhood)
%
%   See also IMDILATE, IMERODE, OFFSETSTREL.

%   Copyright 1993-2015 The MathWorks, Inc.

classdef strel < matlab.mixin.CustomDisplay
    
    %----------------------------------------------------------------------
    properties (Access = protected)

        Type
        nhood
        height
        decomposition
        version = 3;
        
        % precomputed properties
        Flat
        CachedSequence
        
    end
    
    %----------------------------------------------------------------------
    properties (Dependent)
       
        Neighborhood         
        Dimensionality    
        
    end
    
    % public methods ------------------------------------------------------
    methods
        
        %------------------------------------------------------------------
        function se = strel(varargin)
            
            if (nargin == 0)
                % No input arguments --- return empty strel
                %default Type is arbitrary, even for empty strels.
                se.Type = 'arbitrary';

            elseif (nargin == 1) && ...
                (isa(varargin{1}, 'strel') || isa(varargin{1}, 'struct'))
                
                % One strel or struct input --- copy constructor
                % automatically update version by using default                 
                se1 = varargin{1};
                
                se = repmat(strel, [1 numel(se1)]);
                for j = 1:numel(se1)
                    se(j) = strel;
                    se(j).nhood = se1(j).nhood;
                    se(j).height = se1(j).height;
                    se(j).decomposition = se1(j).decomposition;
                    if isequal(se1(j).version, 3)
                        se(j).Type = se1(j).Type;
                    else
                        se(j).Type = 'arbitrary';
                    end                    
                end
                
            else
                [type,params] = ParseInputs(varargin{:});
                
                switch type
                    case 'arbitrary'
                        se = MakeArbitraryStrel(params{:});
                    case 'square'
                        se = MakeSquareStrel(params{:});
                    case 'diamond'
                        se = MakeDiamondStrel(params{:});
                    case 'rectangle'
                        se = MakeRectangleStrel(params{:});
                    case 'octagon'
                        se = MakeOctagonStrel(params{:});
                    case 'line'
                        se = MakeLineStrel(params{:});
                    case 'pair'
                        se = MakePairStrel(params{:});
                    case 'periodicline'
                        se = MakePeriodicLineStrel(params{:});
                    case 'disk'
                        se = MakeDiskStrel(params{:});
                    case 'ball'
                        se = MakeBallStrel(params{:});
                    case 'sphere'
                        se = MakeSphereStrel(params{:});
                    case 'cube'
                        se = MakeCubeStrel(params{:});
                    case 'cuboid'
                        se = MakeCuboidStrel(params{:});
                    otherwise
                        error(message('images:strel:unknownStrelType'));
                end
                se.Type = type;
            end
            
            % precompute Flat
            for j = 1:numel(se)
                se(j).Flat = ~any(se(j).height(:)); %#ok<AGROW>
                se(j).CachedSequence = se(j).decomposition(:);
            end
            
        end

        %------------------------------------------------------------------
        function TF = isequal(se1,se2)
        
            TF = isequal(class(se1), class(se2));
            TF = TF && all(se1.eq(se2));
            
        end
        
        %------------------------------------------------------------------
        function TF = eq(se1,se2)
            
            sizesMatch = isequal(size(se1),size(se2));
            if (~sizesMatch)
                TF = false;
                return
            end

            if isempty(se1) % se2 must be empty too
                TF = true;
                return
            end
            
            % Compare all properties except CachedSequence
            TF = false(1, numel(se1));
            for k = 1:numel(se1)
                TF(k) = ...
                    isequal(se1(k).nhood, se2(k).nhood) && ...
                    isequal(se1(k).height, se2(k).height) && ...
                    isequal(se1(k).decomposition, se2(k).decomposition) && ...
                    isequal(se1(k).version, se2(k).version) && ...
                    isequal(se1(k).Flat, se2(k).Flat);                
            end
            
        end
        
        %------------------------------------------------------------------
        function se2 = reflect(se1)
            %REFLECT Reflect structuring element.
            %   SE2 = REFLECT(SE) reflects a structuring element through its center.
            %   The effect is the same as if you rotated the structuring element's
            %   domain 180 degrees around its center (for a 2-D structuring element).
            %   If SE is an array of structuring element objects, then REFLECT(SE)
            %   reflects each element of SE, and SE2 has the same size as SE.
            %
            %   Example
            %   -------
            %       se = strel([0 0 1; 0 0 0; 0 0 0])
            %       se2 = reflect(se)
            %
            %   See also STREL.
            
            % Testing notes
            % se1:          STREL array; can be empty
            %
            % se2:          STREL array; same size as se1.  Each individual strel
            %               in se2 is the reflection of the corresponding strel
            %               in se1.
            %
            % Note that the reflection operation forces the size of the strel
            % neighborhoods to be odd.  For example:
            % >> se = strel(ones(2,2))
            % se =
            % Flat STREL object containing 4 neighbors.
            %
            % Neighborhood:
            %      1     1
            %      1     1
            % >> se2 = reflect(se)
            % se2 =
            % Flat STREL object containing 4 neighbors.
            %
            % Neighborhood:
            %      1     1     0
            %      1     1     0
            %      0     0     0
            %
            % GETNHOOD(REFLECT(SE)) should be logical.
            
            if length(se1) ~= 1
                % Translate every structuring element in the array.
                se2 = se1;
                for k = 1:numel(se2)
                    se2(k) = reflect(se2(k));
                end
            else
                se2 = strel(se1);
                nhood_local = se2.nhood;
                height_local = se2.height;
                num_dims = ndims(nhood_local);
                subs = cell(1,num_dims);
                size_nhood = size(nhood_local);
                for k = 1:num_dims
                    subs{k} = size_nhood(k):-1:1;
                end
                nhood_local = nhood_local(subs{:});
                height_local = height_local(subs{:});
                new_size = size_nhood + (rem(size_nhood,2) ~= 1);
                if any(new_size > size_nhood)
                    nhood_local = padarray(nhood_local, new_size - size_nhood,0,'post');
                    height_local = padarray(height_local, new_size - size_nhood,0,'post');
                end
                
                se2.nhood = logical(nhood_local);
                se2.height = height_local;
                if ~isempty(se2.decomposition)
                    se2.decomposition = reflect(se1.decompose());
                end
            end
            
        end
        
        %------------------------------------------------------------------
        function se2 = translate(se1,displacement)
            %TRANSLATE Translate structuring element.
            %   SE2 = TRANSLATE(SE,V) translates a structuring element SE in N-D
            %   space.  V is an N-element vector containing the offsets of the
            %   desired translation in each dimension.
            %
            %   Example
            %   -------
            %   Dilating with a translated version of STREL(1) is a way to translate
            %   the input image in space.  This example translates the cameraman.tif
            %   image down and to the right by 25 pixels.
            %
            %       I = imread('cameraman.tif');
            %       se = translate(strel(1), [25 25]);
            %       J = imdilate(I,se);
            %       imshow(I), title('Original')
            %       figure, imshow(J), title('Translated');
            
            % Testing notes
            % se1:          STREL array; can be empty.  Required.
            %
            % v:            Expected to be a vector; if an array is passed in, it
            %               is silently reshaped into row vector.  Must be a double
            %               vector containing integers.  Required.
            %
            % se2:          STREL array; same size as se1.  Each individual strel
            %               in se2 is the translation of the corresponding strel
            %               in se1.
            %
            % Note that the translation operation forces the size of the strel
            % neighborhoods to be odd.  For example:
            % >> se = strel(ones(2,2))
            % se =
            % Flat STREL object containing 4 neighbors.
            %
            % Neighborhood:
            %      1     1
            %      1     1
            % >> se2 = translate(se,[1 1])
            % se2 =
            % Flat STREL object containing 4 neighbors.
            %
            % Neighborhood:
            %      0     0     0     0     0
            %      0     0     0     0     0
            %      0     0     0     0     0
            %      0     0     0     1     1
            %      0     0     0     1     1
            %
            % GETNHOOD(TRANSLATE(SE),V) should be logical.
            %
            % TRANSLATE should work even if length(V) is different than
            % num_dims(getnhood(se)).
            
            if ~isa(se1,'strel')
                error(message('images:translate:wrongType'));
            end
            if ~isa(displacement,'double')
                error(message('images:translate:invalidInputDouble'));
            end
            if any(displacement ~= floor(displacement))
                error(message('images:translate:invalidInput'));
            end
            
            displacement = displacement(:)';
            
            if length(se1) ~= 1
                % Translate every structuring element in the array.
                se2 = se1;
                for k = 1:numel(se2)
                    se2(k) = translate(se2(k), displacement);
                end
            else
                se2 = strel(se1);
                nhood_local = se1.nhood;
                nhood_dims = ndims(nhood_local);
                displacement_dims = length(displacement);
                if (nhood_dims > displacement_dims)
                    displacement = [displacement, zeros(1,nhood_dims - displacement_dims)];
                    num_dims = nhood_dims;
                    size_nhood = size(nhood_local);
                    
                else
                    num_dims = displacement_dims;
                    size_nhood = [size(nhood_local), ones(1,displacement_dims - nhood_dims)];
                end
                
                height_local = se1.height;
                idx = find(nhood_local);
                idx = idx(:);
                sub = cell(1,num_dims);
                [sub{:}] = ind2sub(size_nhood, idx);
                center = floor((size_nhood + 1)/2);
                subs = [sub{:}];
                subs = bsxfun(@minus, subs, center);
                subs = bsxfun(@plus, subs, displacement);              
                max_abs_subs = max(abs(subs),[],1);
                new_size = 2*abs(max_abs_subs) + 1;
                new_center = floor((new_size + 1)/2);
                subs = bsxfun(@plus, subs, new_center);
                for k = 1:num_dims
                    sub{k} = subs(:,k);
                end
                new_idx = sub2ind(new_size, sub{:});
                new_nhood = zeros(new_size);
                new_height = zeros(new_size);
                new_nhood(new_idx) = 1;
                new_nhood = logical(new_nhood);
                new_height(new_idx) = height_local(idx);
                
                se2.nhood = logical(new_nhood);
                se2.height = new_height;
                if (~isempty(se2.decomposition))
                    se2.decomposition(1) = translate(se2.decomposition(1),displacement);
                end
                
            end                     
            
        end
    
        %------------------------------------------------------------------
        function se2 = decompose(obj)
            %decompose Extract sequence of decomposed structuring elements.
            %   SEQ = decompose(SE), where SE is a structuring element array,
            %   returns another structuring element array SEQ containing the
            %   individual structuring elements that form the decomposition of SE.
            %   SEQ is equivalent to SE, but the elements of SEQ have no
            %   decomposition.
            %
            %   Example
            %   -------
            %   STREL uses decomposition for square structuring elements larger than
            %   3-by-3.  Use decompose to extract the decomposed structuring
            %   elements:
            %
            %       se = strel('square',5)
            %       seq = decompose(se)
            %
            %   Use IMDILATE with the 'full' option to see that dilating sequentially
            %   with the decomposed structuring elements really does form a 5-by-5
            %   square:
            %
            %       imdilate(1,seq,'full')
            
            % Testing notes
            % se:          STREL array; individual elements may or may not be
            %              decomposed.
            %
            % seq:         STREL array; individual elements may not be
            % decomposed.
            %              That is, length(decompose(seq(k))) must be 1.  seq
            %              should be a column vector.
                
            if isscalar(obj)
                if ~isempty(obj.CachedSequence)
                    se2 = obj.CachedSequence;
                else
                    if isempty(obj.decomposition)
                        se2 = obj;
                    else
                        % Nested decomposition structures
                        % are not allowed.
                        se2 = obj.decomposition;
                    end
                    
                    obj.CachedSequence = se2;
                end
                
            elseif isempty(obj)
                % A bit of a hack here to return a 1-by-0 strel array.
                se2 = strel;
                se2(1) = [];
                
            else %array of strels.
                
                obj = obj(:);
                se2 = [];
                for k = 1:length(obj)
                    if isempty(obj(k).decomposition)
                        se2 = [se2 obj(k)];
                    else
                        se2 = [se2 obj(k).decomposition];
                    end
                end
                se2 = se2(:);
            end 
        end
        
    end % public methods
    
    % get methods ---------------------------------------------------------
    methods 
        
        %------------------------------------------------------------------
        function dims = get.Dimensionality(obj)
           
            dims = ndims(obj.nhood);
            
        end
        
        %------------------------------------------------------------------
        function Neighborhood = get.Neighborhood(obj)
           
            Neighborhood = obj.nhood;
            
        end 
        
    end
    
    % hidden public methods -----------------------------------------------
    methods (Hidden = true)
        
        %------------------------------------------------------------------
        function height = getheight(se)
            %GETHEIGHT Get height of structuring element.
            %   H = GETHEIGHT(SE) returns an array the same size as GETNHOOD(SE)
            %   containing the height associated with each of the structuring element
            %   neighbors.  H is all zeros for a flat structuring element.
            %
            %   Example
            %   -------
            %       se = strel(ones(3,3),magic(3));
            %       getheight(se)
            
            % Testing notes
            % Syntaxes
            % --------
            % HEIGHT = GETHEIGHT(SE)
            %
            % se:      a 1-by-1 strel array; it may have no neighbors.
            % height:  a double array with the same size as NHOOD = GETNHOOD(SE).
            
            if length(se) ~= 1
                error(message('images:getheight:wrongType'));
            end
            
            height = se.height;
            
        end
        
        %------------------------------------------------------------------
        function [offsets,heights] = getneighbors(se)
            %GETNEIGHBORS Get structuring element neighbor locations and heights.
            %   [OFFSETS,HEIGHTS] = GETNEIGHBORS(SE) returns the relative locations
            %   and corresponding heights for each of the neighbors in the
            %   structuring element object SE.  OFFSETS is a P-by-N array where P is
            %   the number of neighbors in the structuring element and N is the
            %   dimensionality of the structuring element.  Each row of OFFSETS
            %   contains the location of the corresponding neighbor, relative to the
            %   center of the structuring element.  HEIGHTS is a P-element column
            %   vector containing the height of each structuring element neighbor.
            %
            %   Example
            %   -------
            %       se = strel([1 0 1],[5 0 -5])
            %       [offsets,heights] = getneighbors(se)
            
            % Testing notes
            % Syntaxes
            % --------
            % [OFFSETS,HEIGHTS] = GETNEIGHBORS(SE)
            %
            % se:      a 1-by-1 strel array
            %
            % offsets: a num_neighbors-by-num_dims array containing the relative
            %          offsets of each neighbor relative to the center of the
            %          neighborhood.
            %
            % heights: a num_neighbors-by-1 column vector of containing the height
            %          corresponding to each neighbor.
            
            if length(se) ~= 1
                error(message('images:getneighbors:wrongType'));
            end
            
            num_dims = ndims(se.nhood);
            idx = find(se.nhood);
            heights = se.height(idx);
            size_nhood = size(se.nhood);
            center = floor((size_nhood+1)/2);
            subs = cell(1,num_dims);
            [subs{:}] = ind2sub(size_nhood,idx);
            offsets = [subs{:}];
            offsets = reshape(offsets,length(idx),num_dims);
            offsets = bsxfun(@minus, offsets, center);
            
        end
        
        %------------------------------------------------------------------
        function nhood = getnhood(se)
            %GETNHOOD Get structuring element neighborhood.
            %   NHOOD = GETNHOOD(SE) returns the neighborhood associated with the
            %   structuring element SE.
            %
            %   Example
            %   -------
            %
            %       se = strel(eye(5));
            %       nhood = getnhood(se)
            
            % Testing notes
            % Syntaxes
            % --------
            % NHOOD = GETNHOOD(SE)
            %
            % se:       1-by-1 STREL array
            %
            % nhood:    Double array containing 0s and 1s.  Should be logical.
            
            if length(se) ~= 1
                error(message('images:getnhood:wrongType'));
            end
            
            nhood = se.nhood;
            
        end
        
        %------------------------------------------------------------------
        function seq = getsequence(se)          
            seq = se.decompose();
        end        
        
        %------------------------------------------------------------------
        function tf = isflat(se)
            %ISFLAT Return true for flat structuring element.
            %   ISFLAT(SE) returns true (1) if the structuring element SE is flat;
            %   otherwise it returns false (0).  If SE is a STREL array, then TF is
            %   the same size as SE.
            
            % Testing notes            
            % se:           STREL array; can be empty
            %
            % tf:           double logical array, same size as se, containing 0s and
            %               1s.

            tf = false(size(se));
            for k = 1:numel(se)
                tf(k) = se(k).Flat;
            end

        end
        
        %------------------------------------------------------------------
        function se2 = reflect2dFlatStrel(se1)
            % REFLECT2DFLATSTREL Reflect a 2-d, flat structuring element.
            %
            %   SE2 = REFLECT2DFLATSTREL(SE1) reflects a flat,2-D 
            %   structuring element through its center. SE1 may be a STREL 
            %   object or an array of STREL objects. If an array of STREL's
            %   is passed, each STREL should be non-decomposable.
            %
            %   NOTE:
            %   This is a light-weight reflection method that is intended
            %   to be called by @gpuArray/private/morphop.m. Call this
            %   method after strel.decompose().
            
            se2 = strel(se1);
            for k = 1 : numel(se2)
                % Extract the neighborhood and rotate by 180.
                nhood_local = se2(k).nhood;
                size_nhood = size(nhood_local);
                nhood_local = nhood_local(end:-1:1,end:-1:1);
                
                new_size = size_nhood + (rem(size_nhood,2) ~= 1);
                if any(new_size > size_nhood)
                    nhood_local = padarray(nhood_local, new_size - size_nhood,0,'post');
                end
                se2(k).nhood = logical(nhood_local);
                se2(k).height = zeros(new_size);
            end
        end
        
        %------------------------------------------------------------------
        function tf = isdecompositionorthogonal(se)
            
            seq = se.decompose();
            % Algorithm:
            % * Create a matrix P where each row contains the size vector for the nhood of
            %   one strel in the sequence se. (se is already a sequence so we don't need to
            %   call decompose(se) again here.)
            %
            % * Trailing singleton dimensions are added as needed for strels with fewer
            %   dimensions.
            %
            % * Find all values of P~=1. If there's only one of these in each column,
            %   the set of strels se is orthogonal.
            
            num_strels = numel(seq);
            
            P = ones(num_strels,2);
            
            for sInd = 1:num_strels
                nhood_size = size(getnhood(seq(sInd)));
                P(sInd,1:numel(nhood_size)) = nhood_size;
            end
            
            % Fill in trailing singleton dimensions as needed
            P(P==0) = 1;
            
            tf = any( sum(P~=1,1) == 1);

        end
        
        %------------------------------------------------------------------
        function [pad_ul, pad_lr] = getpadsize(se)
            % Required padding 
            
            % Find the array offsets and heights for each structuring element
            % in the sequence.
            seq = decompose(se);
            num_strels = numel(seq);
            offsets    = cell(1,num_strels);
            for sInd = 1:num_strels
                offsets{sInd} = getneighbors(seq(sInd));
            end
            
            if isempty(offsets)
                pad_ul = zeros(1,2);
                pad_lr = zeros(1,2);
                
            else
                num_dims = size(offsets{1},2);
                for k = 2:length(offsets)
                    num_dims = max(num_dims, size(offsets{k},2));
                end
                for k = 1:length(offsets)
                    offsets{k} = [offsets{k} zeros(size(offsets{k},1),...
                        num_dims - size(offsets{k},2))];
                end
                
                pad_ul = zeros(1,num_dims);
                pad_lr = zeros(1,num_dims);
                
                for k = 1:length(offsets)
                    offsets_k = offsets{k};
                    if ~isempty(offsets_k)
                        pad_ul = pad_ul + max(0, -min(offsets_k,[],1));
                        pad_lr = pad_lr + max(0, max(offsets_k,[],1));
                    end
                end
                
            end
        end
        
    end % hidden methods
    
    % display methods -----------------------------------------------------
    methods (Access = protected)
       
        function header = getHeader(obj)
          
            headerStr = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
            
            if isequal(numel(obj), 1)
                
                header = sprintf('%s\n',[headerStr, ' ', getString(message('images:strel:flatObjectHeaderInfo',obj.Type))]);
            else
                header = getHeader@matlab.mixin.CustomDisplay(obj);           
            end

        end        

        function group = getPropertyGroups(~)
            plist = {'Neighborhood','Dimensionality'};
            group = matlab.mixin.util.PropertyGroup(plist,'');
        end
        
    end % display methods

    methods (Static = true)

        %------------------------------------------------------------------
        function b = loadobj(a)

            b = strel(a);           
            
        end 
        
    end

   % ----------------------------------------------------------------------
   methods(Access=private, Static)
       
       %------------------------------------------------------------------
       function name = matlabCodegenRedirect(~)
         name = 'images.internal.coder.strel.StructuringElementHelper';
       end
       
   end
   
end % classdef

%%%
%%% MakeArbitraryStrel
%%%
function se = MakeArbitraryStrel(nhood,height)

se = strel;
se.nhood = nhood ~= 0;
se.height = height;

if (~isempty(nhood) && all(nhood(:)) && ~any(height(:)))
    % Strel is flat with an all-ones neighborhood.  Decide whether to decompose
    % it.
    size_nhood = size(nhood);
    % Heuristic --- if theoretical computation advantage is
    % at least a factor of two, then assume that the advantage
    % is worth the overhead cost of performing dilation or erosion twice.
    advantage = prod(size_nhood) / sum(size_nhood);
    if (advantage >= 2)
        num_dims = ndims(nhood);
        se.decomposition = strel;
        for k = 1:ndims(nhood)
            size_k = ones(1,num_dims);
            size_k(k) = size(nhood,k);
            se.decomposition(k) = strel(ones(size_k));
        end       
    end
end

end

%%%
%%% MakeSquareStrel
%%%
function se = MakeSquareStrel(M)

se = strel(true(M,M));

end

%%%
%%% MakeRectangleStrel
%%%
function se = MakeRectangleStrel(MN)

se = strel(true(MN(1), MN(2)));

end

%%%
%%% MakeDiamondStrel
%%%
function se = MakeDiamondStrel(M)

se = strel;
[rr,cc] = meshgrid(-M:M);
se.nhood = (abs(rr) + abs(cc)) <= M;
se.height = zeros(size(se.nhood));

% Heuristic --- if M > 2, assume computational advantage of decomposition
% is worth the cost of performing multiple dilations (or erosions).
if (M > 2)
    % Compute the logarithmic decomposition of the strel using the method in
    % Rein van den Boomgard and Richard van Balen, "Methods for Fast
    % Morphological Image Transforms Using Bitmapped Binary Images," CVGIP:
    % Models and Image Processing, vol. 54, no. 3, May 1992, pp. 252-254.
    
    n = floor(log2(M));
    se.decomposition = strel([0 1 0; 1 1 1; 0 1 0]);
    for k = 0:(n-1)
        P = 2^(k+1) + 1;
        middle = (P+1)/2;
        nhood = zeros(P,P);
        nhood(1,middle) = 1;
        nhood(P,middle) = 1;
        nhood(middle,1) = 1;
        nhood(middle,P) = 1;
        se.decomposition(end+1) = strel(nhood);
    end
    q = M - 2^n;
    if (q > 0)
        P = 2*q+1;
        middle = (P+1)/2;
        nhood = zeros(P,P);
        nhood(1,middle) = 1;
        nhood(P,middle) = 1;
        nhood(middle,1) = 1;
        nhood(middle,P) = 1;
        se.decomposition(end+1) = strel(nhood);
    end
end

end

%%%
%%% MakeOctagonStrel
%%%
function se = MakeOctagonStrel(M)

% The ParseInputs routine checks to make sure M is a multiple of 3.
k = M/3;
se = strel;

[rr,cc] = meshgrid(-M:M);
se.nhood = abs(rr) + abs(cc) <= M + k;
se.height = zeros(size(se.nhood));

% Compute the decomposition.  To decompose an octagonal strel for M=3k,
% first the strel is decomposed into k strels that each have M=3.  Then,
% each M=3 strel is further (recursively) decomposed into 4 line-segment
% strels.
    
if (k >= 1)
    % Decompose into 4*k strels, each of which have a 3x3 neighborhood.
    a = [0 0 0; 1 1 1; 0 0 0];
    b = a';
    c = eye(3);
    d = rot90(c);
    se.decomposition = strel(a);
    se.decomposition(2) = strel(b);
    se.decomposition(3) = strel(c);
    se.decomposition(4) = strel(d);
    se.decomposition = repmat(se.decomposition, 1, k);
end

end

%%%
%%% MakePairStrel
%%%
function se = MakePairStrel(MN)

se = strel;
size_nhood = abs(MN) * 2 + 1;
se.nhood = false(size_nhood);
center = floor((size_nhood + 1)/2);
se.nhood(center(1),center(2)) = 1;
se.nhood(center(1) + MN(1), center(2) + MN(2)) = 1;
se.height = zeros(size_nhood);

end

%%%
%%% MakeLineStrel
%%%
function se = MakeLineStrel(len,theta_d)

se = strel;

if (len >= 1)
    % The line is constructed so that it is always symmetric with respect
    % to the origin.
    theta = theta_d * pi / 180;
    x = round((len-1)/2 * cos(theta));
    y = -round((len-1)/2 * sin(theta));
    [c,r] = iptui.intline(-x,x,-y,y);
    M = 2*max(abs(r)) + 1;
    N = 2*max(abs(c)) + 1;
    se.nhood = false(M,N);
    idx = sub2ind([M N], r + max(abs(r)) + 1, c + max(abs(c)) + 1);
    se.nhood(idx) = 1;
    se.height = zeros(M,N);
else
    % Do nothing here, which effectively returns the empty strel.
end

end

%%%
%%% MakePeriodicLineStrel
%%%
function se = MakePeriodicLineStrel(p,v)
se = strel;
v = v(:)';
p = (-p:p)';
pp = repmat(p,1,2);
rc = bsxfun(@times, pp, v);
r = rc(:,1);
c = rc(:,2);
M = 2*max(abs(r)) + 1;
N = 2*max(abs(c)) + 1;
se.nhood = false(M,N);
idx = sub2ind([M N], r + max(abs(r)) + 1, c + max(abs(c)) + 1);
se.nhood(idx) = 1;
se.height = zeros(M,N);

end

%%%
%%% MakeDiskStrel
%%%
function se = MakeDiskStrel(r,n)

if (r < 3)
    % Radius is too small to use decomposition, so force n=0.
    n = 0;
end

se = strel;

if (n == 0)
    % Use simple Euclidean distance formula to find the disk neighborhood.  No
    % decomposition.
    [xx,yy] = meshgrid(-r:r);
    nhood = xx.^2 + yy.^2 <= r^2;
    
else
    % Reference for radial decomposition of disks:  Rolf Adams, "Radial
    % Decomposition of Discs and Spheres," CVGIP:  Graphical Models and
    % Image Processing, vol. 55, no. 5, September 1993, pp. 325-332.
    %
    % The specific decomposition technique used here is radial
    % decomposition using periodic lines.  The reference is:  Ronald
    % Jones and Pierre Soille, "Periodic lines: Definition, cascades, and
    % application to granulometries," Pattern Recognition Letters,
    % vol. 17, 1996, pp. 1057-1063.
    
    % Determine the set of "basis" vectors to be used for the
    % decomposition.  The rows of v will be used as offset vectors for
    % periodic line strels.
    switch n
        case 4
            v = [ 1 0
                1 1
                0 1
                -1 1];
            
        case 6
            v = [ 1 0
                1 2
                2 1
                0 1
                -1 2
                -2 1];
            
        case 8
            v = [ 1 0
                2 1
                1 1
                1 2
                0 1
                -1 2
                -1 1
                -2 1];
            
        otherwise
            % This error should have been caught already in ParseInputs.
            error(message('images:getheight:invalidN'));
    end
    
    % Determine k, which is the desired radial extent of the periodic
    % line strels.  For the origin of this formula, see the second
    % paragraph on page 328 of the Rolf Adams paper.
    theta = pi/(2*n);
    k = 2*r/(cot(theta) + 1/sin(theta));
    
    % For each periodic line strel, determine the repetition parameter,
    % rp.  The use of floor() in the computation means that the resulting
    % strel will be a little small, but we will compensate for this
    % below.
    for q = 1:n
        rp = floor(k / norm(v(q,:)));
        if (q == 1)
            se.decomposition = strel('periodicline', rp, v(q,:));
        else
            se.decomposition(q) = strel('periodicline', rp, v(q,:));
        end
    end
    
    % Now dilate the strels in the decomposition together to see how
    % close we came to the desired disk radius.
    
    nhood = imdilate(1, se.decomposition, 'full');
    nhood = nhood > 0;
    [rd,~] = find(nhood);
    M = size(nhood,1);
    rd = rd - floor((M+1)/2);
    max_horiz_radius = max(rd(:));
    radial_difference = r - max_horiz_radius;
    
    % Now we are going to add additional vertical and horizontal line
    % strels to compensate for the fact that the strel resulting from the
    % above decomposition tends to be smaller than the desired size.
    len = 2*(radial_difference-1) + 1;
    if (len >= 3)
        % Add horizontal and vertical line strels.
        se.decomposition(end+1) = strel('line',len,0);
        se.decomposition(end+1) = strel('line',len,90);
        
        % Update the computed neighborhood to reflect the additional strels in
        % the decomposition.
        nhood = imdilate(nhood, se.decomposition(end-1:end), 'full');
        nhood = nhood > 0;
    end
end

se.nhood = nhood;
se.height = zeros(size(nhood));

end

%%%
%%% MakeBallStrel
%%%
function se = MakeBallStrel(r,h,n)

se = strel;

if (r == 0)
    % Make a unit strel.
    se.nhood = true;
    se.height = h;
    
elseif (n == 0)
    % Use Euclidean distance and ellipsoid formulas to construct strel;
    % no decomposition used.
    [xx,yy] = meshgrid(-r:r);
    se.nhood = xx.^2 + yy.^2 <= r^2;
    se.height = h * sqrt(r^2 - min(r^2,xx.^2 + yy.^2)) / r;
    
else
    % Radial decomposition of a sphere.  Reference is the Rolf Adams
    % paper listed above.
    
    % Height profile for each radial line strel is given by a parametric
    % formula of the form (a, g(a)), where a is a function of beta.  See
    % page 331 of the Rolf Adams paper.  Our strategy for using this
    % function is to create a table of (a,g(a)) values and then
    % interpolate into this table.
    beta = linspace(0,pi,100)';
    a = beta - pi/2 - sin(beta).*cos(beta);
    g_a = sin(beta).^2;
    
    % Length of each line strel.
    L = pi*r/n;
    
    % Compute the end-point coordinates of each line strel.
    theta = pi * (0:(n/2 - 1))' / n;
    xy = round(L/2 * [cos(theta) sin(theta)]);
    xy = [xy ; [-xy(:,2) xy(:,1)]];
    
    for k = 1:n
        % For each line strel, compute the x-y coordinates of the elements
        % of the strel, and also compute the corresponding height.
        x = xy(k,1);
        y = xy(k,2);
        [xx,yy] = iptui.intline(0,x,0,y);
        xx = [xx; -xx(2:end)]; %#ok<AGROW>
        yy = [yy; -yy(2:end)]; %#ok<AGROW>
        dist = sqrt(xx.^2 + yy.^2);
        ap = dist*n/r;
        z = h/n * interp1q(a, g_a, ap);
        
        % We could have nan's at the end-points now; replace them by 0.
        z(isnan(z)) = 0;
        
        % Now form neighborhood and height matrices with which we can call
        % strel.
        xmin = min(xx);
        ymin = min(yy);
        M = -2*ymin + 1;
        N = -2*xmin + 1;
        nhood = zeros(M,N);
        height = zeros(M,N);
        row = yy - ymin + 1;
        col = xx - xmin + 1;
        idx = row + M*(col-1);
        nhood(idx) = 1;
        height(idx) = z;
        if (k == 1)
            se.decomposition = strel(nhood,height);
        else
            se.decomposition(k) = strel(nhood,height);
        end
    end
    
    % Now compute the neighborhood and height of the strel resulting the radial
    % decomposition.
    full_height = imdilate(0,se.decomposition,'full');
    full_nhood = isfinite(full_height);
    se.nhood = full_nhood;
    se.height = full_height;
end

end

%%%
%%% MakeSphereStrel
%%%
function se = MakeSphereStrel(r)

    se = strel;

    [x,y,z] = meshgrid(-r:r,-r:r,-r:r);
    se.nhood =  ( (x/r).^2 + (y/r).^2 + (z/r).^2 ) <= 1;
    se.height = zeros(size(se.nhood));

end

%%%
%%% MakeCubeStrel
%%%
function se = MakeCubeStrel(width)

    se = strel(true(width,width,width));

end

%%%
%%% MakeCuboidStrel
%%%
function se = MakeCuboidStrel(XYZ)

se = strel(true(XYZ(1), XYZ(2), XYZ(3)));

end

%%%
%%% ParseInputs
%%%
function [type,params] = ParseInputs(varargin)

default_ball_n = 8;
default_disk_n = 4;

narginchk(1, 4);

if ~ischar(varargin{1})
    type = 'arbitrary';
    params = varargin;
else
    type = varargin{1};
    params = varargin(2:end);
    
    valid_strings = {'arbitrary'
        'square'
        'diamond'
        'rectangle'
        'octagon'
        'line'
        'disk'
        'sphere'
        'cube'
        'cuboid'};

    if strncmpi('ball', type, numel(type))
        type = 'ball';
    elseif strncmpi('pair', type, numel(type))
        type = 'pair';
    elseif strncmpi('periodicline', type, numel(type))
        type = 'periodicline';
    else
        type = validatestring(varargin{1}, valid_strings, 'strel', ...
            'STREL_TYPE', 1);
    end
end

num_params = numel(params);

switch type
    case 'arbitrary'
        if (num_params < 1)
            error(message('images:strel:tooFewInputs',type))
        end
        if (num_params > 2)
            error(message('images:strel:tooManyInputs',type))
        end
        % Check validity of the NHOOD argument.
        nhood = params{1};
        validateattributes(nhood, {'numeric', 'logical'}, {'real'}, 'strel', ...
            'NHOOD', 2);
        
        % Check validity of the HEIGHT argument.
        if num_params >= 2
            height = params{2};
            validateattributes(height, {'double'}, {'real', 'nonnan'}, 'strel', ...
                'HEIGHT', 3);
            if ~isequal(size(height), size(nhood))
                error(message('images:strel:sizeMismatch'))
            end
        else
            params{2} = zeros(size(nhood));
        end
        
    case 'square'
        if (num_params < 1)
            error(message('images:strel:tooFewInputs',type))
        end
        if (num_params > 1)
            error(message('images:strel:tooManyInputs',type))
        end
        M = params{1};
        validateattributes(M, {'double'}, {'scalar' 'integer' 'real' 'nonnegative'}, ...
            'strel', 'SIZE', 2);
        
    case 'diamond'
        if (num_params < 1)
            error(message('images:strel:tooFewInputs',type))
        end
        if (num_params > 1)
            error(message('images:strel:tooManyInputs',type))
        end
        M = params{1};
        validateattributes(M, {'double'}, {'scalar' 'integer' 'nonnegative'}, ...
            'strel', 'SIZE', 2);
        
    case 'octagon'
        if (num_params < 1)
            error(message('images:strel:tooFewInputs',type))
        end
        if (num_params > 1)
            error(message('images:strel:tooManyInputs',type))
        end
        M = params{1};
        validateattributes(M, {'double'}, {'scalar' 'integer' 'nonnegative'}, ...
            'strel', 'SIZE', 2);
        if rem(M,3) ~= 0
            error(message('images:strel:notMultipleOf3'))
        end
        
    case 'rectangle'
        if (num_params < 1)
            error(message('images:strel:tooFewInputs',type))
        end
        if (num_params > 1)
            error(message('images:strel:tooManyInputs',type))
        end
        MN = params{1};
        validateattributes(MN, {'double'}, {'vector' 'real' 'integer' 'nonnegative'}, ...
            'strel', 'SIZE', 2);
        if numel(MN) ~= 2
            error(message('images:strel:badSizeForRectangle'))
        end
        
    case 'pair'
        if (num_params < 1)
            error(message('images:strel:tooFewInputs',type))
        end
        if (num_params > 1)
            error(message('images:strel:tooManyInputs',type))
        end
        RC = params{1};
        validateattributes(RC, {'double'}, {'vector' 'real' 'integer'}, ...
            'strel', 'OFFSET', 2);
        if numel(RC) ~= 2
            error(message('images:strel:badOffsetsForPair'))
        end
        
    case 'line'
        if (num_params < 2)
            error(message('images:strel:tooFewInputs',type))
        end
        if (num_params > 2)
            error(message('images:strel:tooManyInputs',type))
        end
        len = params{1};
        validateattributes(len, {'double'}, {'scalar' 'real' 'finite' 'nonnegative'}, ...
            'strel', 'LEN', 2);
        
        deg = params{2};
        validateattributes(deg, {'double'}, {'scalar' 'real' 'finite'}, 'strel', 'DEG', 3);
        
    case 'periodicline'
        if (num_params < 2)
            error(message('images:strel:tooFewInputs',type))
        end
        if (num_params > 2)
            error(message('images:strel:tooManyInputs',type))
        end
        p = params{1};
        validateattributes(p, {'double'}, {'scalar' 'real' 'integer' 'nonnegative'}, ...
            'strel', 'P', 2);
        
        v = params{2};
        validateattributes(v, {'double'}, {'vector' 'real' 'integer'}, 'strel', ...
            'V', 3);
        if numel(v) ~= 2
            error(message('images:strel:wrongSizeForV'))
        end
        
    case 'disk'
        if (num_params < 1)
            error(message('images:strel:tooFewInputs',type))
        end
        if (num_params > 2)
            error(message('images:strel:tooManyInputs',type))
        end
        
        r = params{1};
        validateattributes(r,{'double'}, {'scalar' 'real' 'integer' 'nonnegative'}, ...
            'strel', 'R', 2);
        
        if (num_params < 2)
            params{2} = default_disk_n;
        else
            n = params{2};
            validateattributes(n, {'double'}, {'scalar' 'real' 'integer'}, ...
                'strel', 'N', 3);
            if ((n ~= 0) && (n ~= 4) && (n ~= 6) && (n ~= 8))
                error(message('images:strel:invalidN'))
            end
        end
        
    case 'ball'
        if (num_params < 2)
            error(message('images:strel:tooFewInputs',type))
        end
        if (num_params > 3)
            error(message('images:strel:tooManyInputs',type))
        end
        
        r = params{1};
        validateattributes(r, {'double'}, {'scalar' 'real' 'integer' 'nonnegative'}, ...
            'strel', 'R', 2);
        
        h = params{2};
        validateattributes(h, {'double'}, {'scalar' 'real'}, 'strel', 'H', 3);
        
        if (num_params < 3)
            params{3} = default_ball_n;
        else
            n = params{3};
            validateattributes(n, {'double'}, {'scalar' 'integer' 'nonnegative' ...
                'even'}, 'strel', 'N', 4);
        end
        
    case 'sphere'        
        if (num_params < 1)
            error(message('images:strel:tooFewInputs',type))
        end
        if (num_params > 1)
            error(message('images:strel:tooManyInputs',type))
        end
        
        r = params{1};        
        validateattributes(r, {'double'}, {'scalar' 'integer' 'nonnegative'}, ...
            'strel', 'R', 2);
        
    case 'cube'        
        if (num_params < 1)
            error(message('images:strel:tooFewInputs',type))
        end
        if (num_params > 1)
            error(message('images:strel:tooManyInputs',type))
        end
        
        width = params{1};        
        validateattributes(width, {'double'}, {'scalar' 'integer' 'nonnegative'}, ...
            'strel', 'WIDTH', 2);
        
    case 'cuboid'        
        if (num_params < 1)
            error(message('images:strel:tooFewInputs',type))
        end
        if (num_params > 1)
            error(message('images:strel:tooManyInputs',type))
        end

        XYZ = params{1};
        validateattributes(XYZ, {'double'}, {'vector' 'integer' 'nonnegative'}, ...
            'strel', 'SIZE', 2);
        if numel(XYZ) ~= 3
            error(message('images:strel:badSizeForCuboid'))
        end
        
    otherwise
        % This code should be unreachable.
        error(message('images:strel:unrecognizedStrelType'))
end

end

