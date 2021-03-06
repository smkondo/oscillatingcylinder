### Authors: Nathaniel Ruhl and Shinya Kondo
##### Final Project for MEAM 536: Viscous Fluid Flow. 
##### Institution: University of Pennsylvania

This repository contains MATLAB code, animation movies, and figures that are made to analyze the flow behind a transversely oscillating cylinder. ANSYS was used to perform the CFD simulations. CFD data is available [here](https://drive.google.com/drive/folders/1rjoNhZw0ePAcFWGrpETHtA72gybqfsyj?usp=sharing). The scripts readData.m and readVorticityData.m provide functions that read the velocity and vorticity data from these files. The contents of the CFD data folder is as follows:

- StationaryCylinderFinal/: Data for the stationary cylinder
- moving1final/: Data for frequency ratio of 1
- moving2final/: Data for frequency ratio of 1.5
- moving3final/: Data for frequency ratio of 0.5
- Coefficients of drag and lift are given in the ".txt" files

----------------------------------
###### Contents of movies/:
- stationaryfinal.avi: Streamwise velocity contour for stationary cylinder
- moving1final.avi: Streamwise velocity contour for frequency ratio of 1
- moving2final.avi: Streamwise velocity contour for frequency ratio of 1.5
- moving3final.avi: Streamwise velocity contour for frequency ratio of 0.5
- ROM.avi: Reduced-order model for the stationary cyllinder
- ROM2.avi: Reduced-order model for frequency ratio 1
- vorticity3.avi: Vorticity contour for frequency ratio 1.5
- vorticitymoving1.avi and vorticitymoving1slow.avi: Vorticity contour for frequency ratio 1
- vorticitystationary.avi: Vorticity contour for stationary cylinder





