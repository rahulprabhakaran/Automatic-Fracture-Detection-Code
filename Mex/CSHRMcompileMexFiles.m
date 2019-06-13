%Compile necessary MEX files

cd  Mex
mex CSHRMgetEdgesAndTangentOrientationsMex.cpp
mex CSHRMgetRidgesAndTangentOrientationsMex.cpp
mex CSHRMgetCurvatureMex.cpp

%  Copyright (c) 2016. Rafael Reisenhofer
%
%  Part of CoShREM Toolbox v1.1
%  Built Mon, 11/01/2016
%  This is Copyrighted Material