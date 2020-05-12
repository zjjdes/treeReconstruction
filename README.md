# treeReconstruction
This repository contains the data and code to reproduce the experiment in the paper titled **3D Reconstruction of Internal Wood Decay Using Photogrammetry and Sonic Tomography** (Authors: Junjie Zhang, Kourosh Khoshelham).

## Paper Abstract
Knowledge of deteriorations within tree trunks is critical for arborists to conduct individual tree health assessment. Sonic tree tomography, a non-destructive technique using sound waves, has been widely used to estimate the size and shape of the internal decay based on sound wave velocity variations. However, sonic tree tomography has commonly been applied to 2D horizontal or vertical cross sections and its accuracy is questionable due to the poor approximation of the shape of the cross section. This paper proposes an integration of close-range photogrammetry and sonic tomography to enable accurate reconstruction of the exterior and interior of the tree trunk in 3D. The internal wood quality is represented by the spatially interpolated sound wave velocities, using the time of flight of the sound waves and the coordinates of the acoustic sensors obtained from the photogrammetric model. Experimental results show that the proposed approach provides a realistic 3D visualisation of the size, shape and location of the internal deteriorations.

## Data
The following files contain the input data for the experiment used in the paper (data was collected by Junjie Zhang and David Galwey):
- `Runtimes.mat`
- `Markers.mat`
- `trunk.ply`

## Code
The following files contain the code for the experiment used in the paper and should be executed in the indicated order:
1. `control_points.m`
2. `truncate_control_points.m`
3. `semi_variogram.m`
4. `interpolation.m`

while `svfun.m` contains the semi-variance function which should be obtained using `cftool` in a step in `semi_variogram.m`. Some settings can be modified in these scripts to observe changes in the results.
