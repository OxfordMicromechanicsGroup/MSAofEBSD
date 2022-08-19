# MSAofEBSD
multivariate statistical analysis of EBSD patterns using MatLab

Matlab scripts to implement analysis described in detail through paper below. 

MSAofEBSD_HDF5_v01.m
Takes set of EBSD patterns forming a map of sample, and indentifies similar patterns which are clasified into distinct grains.
Classification is through PCA analysis, VARIMAX solution, or k-means clustering (slower).
EBSD data provided in hdf5 format, and user input of number of grains to segment map into.
Characteristic EBSD patterns are generated for each grains

Driver_PatMatch_folder_steel_1.m
The charactersitic patterns can then be indexed by comparing to a template libray generated by interopolating intensities in a simulated master pattern generated by Bruker Dynamics (or other code).  The master pattern is supplied as intensities on a sterographic projection. 
The Mtex plugin is required for this part of the code.




 "Applications of multivariate statistical methods and simulation libraries to analysis of electron backscatter diffraction and transmission Kikuchi diffraction datasets" 
Angus Wilkinson, David Collins, Yevhen Zayachuk, Rajesh Korla, Arantxa Vilalta-Clemente
Ultramicroscopy     -> https://doi.org/10.1016/j.ultramic.2018.09.011
arXiv (open access) -> https://arxiv.org/abs/1806.02087

Abstract:
Multivariate statistical methods are widely used throughout the sciences, including microscopy, however, their utilisation for analysis of electron backscatter diffraction (EBSD) data has not been adequately explored. The basic aim of most EBSD analysis is to segment the spatial domain to reveal and quantify the microstructure, and links this to knowledge of the crystallography (eg crystal phase, orientation) within each segmented region. Two analysis strategies have been explored; principal component analysis (PCA) and k-means clustering. The intensity at individual (binned) pixels on the detector were used as the variables defining the multidimensional space in which each pattern in the map generates a single discrete point. PCA analysis alone did not work well but rotating factors to the VARIMAX solution did. K-means clustering also successfully segmented the data but was computational more expensive. The characteristic patterns produced by either VARIMAX or k-means clustering enhance weak patterns, remove pattern overlap, and allow subtle effects from polarity to be distinguished. Combining multivariate statistical analysis (MSA) approaches with template matching to simulation libraries can significantly reduce computational demand as the number of patterns to be matched is drastically reduced. Both template matching and MSA approaches may augment existing analysis methods but will not replace them in the majority of applications.
