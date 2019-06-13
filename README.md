# Supplementary_Code_SolidEarth 
This repository contains MATLAB scripts that implement a method to extract digitized fractures from images of fractured rocks. The method is described in the manuscript titled ,"An automated fracture trace detection technique using the complex shearlet transform"(se-2019-104) that is submitted to Solid Earth Journal (https://www.solid-earth.net/). 

The MATLAB files in this repository consists of four MATLAB scripts that perform the following functions:

  1. Ridge_Ensemble_Generator.m 
                                -> creates a set of shearlet systems based on user-defined ranges and saves them as *.mat files
                                -> reads a set of images of fractured rocks
                                -> computes a number of ridge realizations based on a user-defined range of ridge extraction parameters
                                   and adds the detected ridges into a single image and saves it. This is referred to as ridge ensemble.
  2. Ridge_Ensemble_Reader.m 
                                -> reads the ridge ensemble image and applies a non-linear, sigmoid function in order to filter away 
                                   possible false positives and retain probable fracture ridges, by means of a user-defined threshold
                                -> saves the highly probable, fracture ridge map as a binarized image                               
  3. Ridge_Post_Processing.m
                                -> reads the binarized ridge image files
                                -> performs a segmentation operation using Otsu thresholding that removes isolated and very small pixel
                                   clusters
                                -> skeletonizes the ridge clusters into thinned pixel reresentations preserving the branch points
                                -> fits polylines to the skeletonized ridge clusters and saves as tables (*.mat files) 
  4. Polyline_to_Shape.m   
                                -> converts the polylines into shape files
                                -> georeferences the polylines (if georeferencing information is available)
                                -> performs rotation, scaling, and translation of polylines (if necessary)
                                -> performs line simplification of polylines using the Douglas-Peucker algorithm (if necessary)
                                -> converts polylines into shapefiles structures and saves as shapefiles

The "Examples" folder contains images which are referred to in the manuscript. These are provided so that the user of this repository
can replicate the results in the manuscript. The images are sourced from the following sources:
   
1. Bingie_Bingie: Thiele, S. T., Grose, L., Samsu, A., Micklethwaite, S., Vollgger, S. A., and Cruden, A. R.: Rapid, semi-automatic      fracture and contact mapping for point clouds, images and geophysical data, Solid Earth, 8 (6), 1241–1253, https://doi.org/10.5194/se-8-1241-2017, 2017.

2. Core_Slice: Dwarkasing, A.: 3D Fracture Analyses of Various Rock Samples through X-Ray Micro-Tomography, Master’s thesis, Delft University of Technology, http://resolver.tudelft.nl/uuid:aefe6746-9788-4c93-ac4c-f80250e6a12c, 2016.

3. Parmelan_Tile: Prabhakaran, R., Bruna, P.-O., Bertotti, G., Smeulders, D., and Meda, M.: Fracture Network Patterns from the Parmelan Anticline, France),4TU Centre for Research Data. Dataset, https://doi.org/10.4121/uuid:3f5e255f-edf7-441f-89f2-1adc7ac2f7d1, 2019b.


The "Dependencies" folder contains MATLAB functions that are sourced from the open-source MATLAB toolboxes, CosHREM and Geom2D, and which are used in the MATLAB scripts.

The "Mex" folder contains C++ executables that are used by the MATLAB scripts. They are from the open-source Complex Shearlet Toolbox (CoSHREM). Please refer to : "Reisenhofer, R., Kiefer, J., and King, E. J.: Shearlet-based detection of flame fronts, Experiments in Fluids, 57, 41, https://doi.org/10.1007/s00348-016-2128-6, 2016" for a detailed description of the CoSHREM toolbox.


