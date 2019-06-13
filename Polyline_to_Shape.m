%-------------------------------------------------------------------------%
% This script converts polyline files that are in *.mat format into       %
% georeferenced shape files (if there is georeferencing information)      %
% which can then be viewed in QGIS / ArcGIS. The polylines can be         %
% simplified before writing as a shapefile using the Douglas-Peucker line %
% simplification algorithm.                                                 
%-------------------------------------------------------------------------%
clc
clear all
close all
format compact

%%  select polyline mat files 
    disp('Select multiple mat files (they must be in a single folder)');
    [InFileListShort, pathname] = uigetfile('*.jpg;*.tif;*.png;*.mat','Select input image files','MultiSelect','on');
    
    if not(iscell(InFileListShort))
       InFileListShort = {InFileListShort};
    end
    
    replicatePath = repmat(cellstr(pathname),size(InFileListShort));
    InFileList = strcat(replicatePath,InFileListShort);    
    InFileList=InFileList';
    
%%  select the world files if there is georeferencing information
    disp('Select world files (they must be in a single folder)');
    [InFileListShort2, pathname2] = uigetfile('*.jpg;*.tfw;*.png;*.mat','Select input image files','MultiSelect','on');
    
    if not(iscell(InFileListShort2))
       InFileListShort2 = {InFileListShort2};
    end
    
    replicatePath2 = repmat(cellstr(pathname2),size(InFileListShort2));
    InFileList2 = strcat(replicatePath2,InFileListShort2);    
    InFileList2=InFileList2';
    
%% creates a list of files which have polylines...some tiles of the original
% orthomosaic (on the boundary) do not have no detected features
% creates a cell array of the fractures detected in each tile

% an aspect ratio can be applied here to adjust the polylines
aspect_ratio = 1/1;
for i=1:length(InFileList)
   load(InFileList{i});   
   [table_size,~] = size(Poly_Points_Table);
    for j=1:table_size
     Poly_Points_Table.Poly_Points{j,1}(:,1)=Poly_Points_Table.Poly_Points{j,1}(:,1).*aspect_ratio;
     Poly_Points_Table.Poly_Points{j,1}(:,2)=Poly_Points_Table.Poly_Points{j,1}(:,2)./aspect_ratio;
     i  
         
    end
   PolyLines{i,1} = table2cell(Poly_Points_Table);
   
end

%% creating an N1N2 matrix
m=length(PolyLines{1,1});

% creating row indices matrix
N1N2=zeros(length(PolyLines),2);
N1N2(1,1)=1;
N1N2(1,2)=length(PolyLines{1,1});
for i=2:length(PolyLines)
  N1N2(i,1)=  N1N2(i-1,2)+1;
  N1N2(i,2)=  length(PolyLines{i,1})+N1N2(i-1,2);   
end 

%% Georeferencing the PolyLines
PolyLines_Georeferenced = PolyLines;

for i=1:length(PolyLines)
   
    R = worldfileread(InFileList2{i}, 'planar', [938 1062]);
    
 % calculating the corner points of the orthotile (in lat-long )
 Lr_x = R.XWorldLimits(1);
 Lr_y = R.YWorldLimits(1);
 
 % converting the pixel values into lat-long and then setting the origin
 % based on the corner values
 for j=1:length(PolyLines{i,1})
   PolyLines_Georeferenced{i,1}{j,1}(:,1)= PolyLines_Georeferenced{i,1}{j,1}(:,1).*...
       R.CellExtentInWorldX + Lr_x;
   
   
   PolyLines_Georeferenced{i,1}{j,1}(:,2)= PolyLines_Georeferenced{i,1}{j,1}(:,2).*...;   
       R.CellExtentInWorldY + Lr_y;
     
 end  
 clearvars  Lr_x Lr_y R
 i
end




%% Rotating the Georeferenced PolyLines
PolyLines_Georeferenced_Rotated = PolyLines_Georeferenced;
j=1;
for i=1:length(PolyLines_Georeferenced) 
    %[~,R] = geotiffread(InFileList2{i});  
    R = worldfileread(InFileList2{i}, 'planar', [938 1062]);
    for k=1:length(PolyLines_Georeferenced{j,1})  
%     C=PolyLines_Georeferenced{j,1};  

      % the following step flips the shapefile about the x-axis
        %PolyLines_Georeferenced_Rotated{j,1}{k,1}(:,2) = PolyLines_Georeferenced{j,1}{k,1}(:,2)*-1 + R.YWorldLimits(1,1) + R.YWorldLimits(1,2) ; 
    
      % the following step flips the shapefile about the y-axis
        PolyLines_Georeferenced_Rotated{j,1}{k,1}(:,1) = PolyLines_Georeferenced{j,1}{k,1}(:,1)*-1 + R.XWorldLimits(1,1) + R.XWorldLimits(1,2) ; 
    end
     j=j+1;          
     disp(i)
end



%% Rotating non-Georeferenced Polylines
PolyLines_Rotated = PolyLines;
j=1;
for i=1:length(PolyLines) 
    
    for k=1:length(PolyLines{j,1})  
%     C=PolyLines_Georeferenced{j,1};  

      % rotating about the x-axis
      %PolyLines_Rotated{j,1}{k,1}(:,2) = PolyLines{j,1}{k,1}(:,2)*-1 + 1000 + 0;
      
      % rotating about the y-axis
       PolyLines_Rotated{j,1}{k,1}(:,1) = PolyLines{j,1}{k,1}(:,1)*-1 + 0 + 1000;
      
      % translating down the y-axis
      %PolyLines_Rotated{j,1}{k,1}(:,2) = PolyLines{j,1}{k,1}(:,2) - 1000;
    end  
     j=j+1;          
     disp(i)
end


%% Creating the ShapeFile Structure (not georeferenced)
j=1;
%P = PolyLines_Rotated;
P = PolyLines;

for i=1:length(P) 
  
  C=P{i,1};
  if i>1
    m=length(C)+length(P{i-1,1});
  end
  z=1;
  for k=N1N2(i,1):N1N2(i,2)  
    [PolyLines_Shape(k).Tile_ID] = i;
%     [PolyLines_Shape(k).Ortho_ID] = P{i,2};
%      [PolyLines_Shape(k).BoundingBox] = [R{P{i,2},1}.LongitudeLimits(1,1) ...
%                                         R{P{i,2},1}.LatitudeLimits(1,1); ...
%                                         R{P{i,2},1}.LongitudeLimits(1,2) ...
%                                         R{P{i,2},1}.LatitudeLimits(1,2)];
    [PolyLines_Shape(k).Geometry] = 'PolyLine';
    [PolyLines_Shape(k).Polyline_ID] = k;
    [PolyLines_Shape(k).Polyline_in_Ortho_ID] = z;
    % if C is a cell array
    % [PolyLines_Shape(k).X] =  C{k-j+1,1}(:,1);
    % if C is not a cell array
    [PolyLines_Shape(k).X] =  C(k,1);
    % if C is a cell array
    % [PolyLines_Shape(k).Y] =  C{k-j+1,1}(:,2);   
    % if C is not a cell array
    [PolyLines_Shape(k).Y] =  C(k,2); 
    z=z+1;
  end 
  j=length(PolyLines_Shape)+1;
  
  clearvars C; clearvars z;
end 

%% Creating the Shapefile Structure (not georeferenced, but rotated)
j=1;
%P = PolyLines;
%P = PolyLines_Rotated;
P = P_Simplified;

for i=1:length(P) 
  
  C=P{i,1};
  if i>1
    m=length(C)+length(P{i-1,1});
  end
  z=1;
  for k=N1N2(i,1):N1N2(i,2)  
    [PolyLines_Shape(k).Tile_ID] = i;
%     [PolyLines_Shape(k).Ortho_ID] = PolyLines{i,2};
%      [PolyLines_Shape(k).BoundingBox] = [R{PolyLines_Georeferenced_Rotated{i,2},1}.LongitudeLimits(1,1) ...
%                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LatitudeLimits(1,1); ...
%                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LongitudeLimits(1,2) ...
%                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LatitudeLimits(1,2)];
    [PolyLines_Shape(k).Geometry] = 'PolyLine';
    [PolyLines_Shape(k).Polyline_ID] = k;
    [PolyLines_Shape(k).Polyline_in_Ortho_ID] = z;
    % if C is a cell array
     [PolyLines_Shape(k).X] =  C{k-j+1,1}(:,1);
    % if C is not a cell array
    %[PolyLines_Shape(k).X] =  C(k,1);
    % if C is a cell array
     [PolyLines_Shape(k).Y] =  C{k-j+1,1}(:,2);   
    % if C is not a cell array
    %[PolyLines_Shape(k).Y] =  C(k,2); 
    z=z+1;
  end 
  j=length(PolyLines_Shape)+1;
  
  clearvars C; clearvars z;
end 



%% Creating the ShapeFile Structure (georeferenced and rotated)
j=1;
%P = PolyLines_Georeferenced_Rotated;
%P = PolyLines_Georeferenced;
P = P_Simplified;
for i=1:length(P) 
  
  C=P{i,1};
  if i>1
    m=length(C)+length(P{i-1,1});
  end
  z=1;
  for k=N1N2(i,1):N1N2(i,2)  
    [PolyLines_Shape(k).Tile_ID] = i;
%     [PolyLines_Shape(k).Ortho_ID] = PolyLines{i,2};
%      [PolyLines_Shape(k).BoundingBox] = [R{PolyLines_Georeferenced_Rotated{i,2},1}.LongitudeLimits(1,1) ...
%                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LatitudeLimits(1,1); ...
%                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LongitudeLimits(1,2) ...
%                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LatitudeLimits(1,2)];
    [PolyLines_Shape(k).Geometry] = 'PolyLine';
    [PolyLines_Shape(k).Polyline_ID] = k;
    [PolyLines_Shape(k).Polyline_in_Ortho_ID] = z;
    % if C is a cell array
     [PolyLines_Shape(k).X] =  C{k-j+1,1}(:,1);
    % if C is not a cell array
    %[PolyLines_Shape(k).X] =  C(k,1);
    % if C is a cell array
    [PolyLines_Shape(k).Y] =  C{k-j+1,1}(:,2);   
    % if C is not a cell array
%     [PolyLines_Shape(k).Y] =  C(k,2); 
    z=z+1;
  end 
  j=length(PolyLines_Shape)+1;
  disp(i)
  clearvars C; clearvars z;
end 

%% Creating ShapeFile Structures and writing shapefile for each image separately

P = P_Simplified;
outfolder='D:\Github_Test\';

% j=1;
for i=1:length(P)
  tic  
  filename=InFileListShort{i};
  filename=filename(1:length(filename)-4);
  outfilename=strcat(filename,'.shp');

  C=P{i,1};
%   if i>1
%     m=length(C)+length(PolyLines_Georeferenced_Rotated{i-1,1});
%   end
  z=1;
  for k=1:length(C)  
    [PolyLines_Shape(k).Tile_ID] = i;
%     [PolyLines_Shape(k).Ortho_ID] = PolyLines{i,2};
%      [PolyLines_Shape(k).BoundingBox] = [R{PolyLines_Georeferenced_Rotated{i,2},1}.LongitudeLimits(1,1) ...
%                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LatitudeLimits(1,1); ...
%                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LongitudeLimits(1,2) ...
%                                         R{PolyLines_Georeferenced_Rotated{i,2},1}.LatitudeLimits(1,2)];
    [PolyLines_Shape(k).Geometry] = 'PolyLine';
    [PolyLines_Shape(k).Polyline_ID] = k;
    [PolyLines_Shape(k).Polyline_in_Ortho_ID] = z;
    % if C is a cell array
     %[PolyLines_Shape(k).X] =  C{k-j+1,1}(:,1);
    % if C is not a cell array
    [PolyLines_Shape(k).X] =  C{k,1}(:,1);
    % if C is a cell array
    %[PolyLines_Shape(k).Y] =  C{k-j+1,1}(:,2);   
     %if C is not a cell array
     [PolyLines_Shape(k).Y] =  C{k,1}(:,2); 
    z=z+1;
  end 
  j=length(PolyLines_Shape)+1;
  shapewrite(PolyLines_Shape,strcat(outfolder,outfilename));
  disp(i)
  clearvars C  PolyLines_Shape;
  toc
end 


%% performing line simplification on the fitted curves using Douglas-Pecker
% algorithm
%epsilon = 1E-10;
epsilon = 1E-02;

%P_Simplified = PolyLines_Georeferenced_Rotated;
%P_Simplified = PolyLines_Georeferenced;
%P_Simplified = PolyLines_Rotated;
P_Simplified = PolyLines_Rotated;
for i=1:length(P_Simplified)
    tic
 for j=1:length(P_Simplified{i,1})
    Points = P_Simplified{i,1}{j,1}';
    Simplified_Points = DouglasPeucker(Points,epsilon);
    P_Simplified{i,1}{j,1} = Simplified_Points';
 end
 disp(i)
 toc
end



%% Writing shapefile
tic
%shapewrite(PolyLines_Shape,'D:\PhD\Automatic_Detection\Core_Fractures\Shape_Files\xz769_Rotated_Simplified.shp');
shapewrite(PolyLines_Shape,'D:\PhD\Automatic_Detection\Probabilistic_Edges\Ortho_1_107\Ortho_107_New\Shapefiles_Georef\Tile_107_Georef.shp');
toc

%% Classifying Polylines by length
for i=1:length(PolyLines{1,1})
  polyline_length = 0;  
  for j=1:100 
    polyline_length = polyline_length + Lengths2D([PolyLines{1,1}{i,1}(j,1) PolyLines{1,1}{i,1}(j,2) PolyLines{1,1}{i,1}(j+1,1) PolyLines{1,1}{i,1}(j+1,2)]);   
  end  
  Poly_Length(i,1) = polyline_length;
end
%%
for i=1:length(New_PolyLines)
  polyline_length = 0;  
  for j=1:100 
    polyline_length = polyline_length + Lengths2D([New_PolyLines{i,1}(j,1) New_PolyLines{i,1}(j,2) New_PolyLines{1,1}(j+1,1) New_PolyLines{i,1}(j+1,2)]);   
  end  
  Poly_Length(i,1) = polyline_length;
end


%% Classifying Polylines by Azimuth
for i=1:length(PolyLines{1,1})   
Poly_Angle(i,1) = Angles2D([PolyLines{1,1}{i,1}(1,1) PolyLines{1,1}{i,1}(1,2) PolyLines{1,1}{i,1}(101,1) PolyLines{1,1}{i,1}(101,2)]);   
  if Poly_Angle(i,1)<0
    Poly_Angle(i,1)=Poly_Angle(i,1)+180;
  end    
end
%% 
% creating a new polylines cell array based on length threshold
clearvars New_Polylines_2
j=1;
for i=1:length(New_PolyLines)
  if Poly_Length(i,1)>5000
      New_PolyLines_2{j,1} =New_PolyLines{i,1};
      j=j+1;
  end          
end

%% creating a new polylines cell array based on angle
j=1;
for i=1:length(PolyLines{1,1})
  if (74<Poly_Angle(i,1)) && (Poly_Angle(i,1)>73)
      New_PolyLines{j,1} =PolyLines{1,1}{i,1};
      j=j+1;
  end          
end

%%  Displaying the classified polylines i.e. New_PolyLines
figure(1)
imshow(cdata); hold on;
for i=1:length(New_PolyLines)
    scatter(New_PolyLines{i,1}(:,1),New_PolyLines{i,1}(:,2),'b.')
    grid on
    pbaspect([1 1 1])
    hold on   
end    

%% Displaying the classified polylines i.e. New_PolyLines_2
figure(1)
imshow(cdata); hold on;
for i=1:length(New_PolyLines_2)
    scatter(New_PolyLines_2{i,1}(:,1),New_PolyLines_2{i,1}(:,2),'b.')
    grid on
    pbaspect([1 1 1])
    hold on   
end 

%% Displaying all polylines 
figure(1)
imshow(cdata); hold on;
for i=1:length(PolyLines{1,1})
    scatter(PolyLines{1,1}{i,1}(:,1),PolyLines{1,1}{i,1}(:,2),'y.')
    grid on
    pbaspect([1 1 1])
    hold on   
end  
