# Simulate confocal fluorescence microscopy and evaluate CaCLEAN performance

This repository contains scripts is used to simulate confocal fluorescence microscopy images
from the results of a [finite element reaction-diffusion model](https://github.com/dladd/CardiacCalcium_FiniteElement) of calcium
diffusing in the intracellular space of a cardiomyocyte. The simulated images
are then used to evaluate the performance of the CaCLEAN algorithm, which seeks
to identify the locations of calcium release sites in experimental imaging of
live cardiomyocytes. 

For more information, please see the corresponding manuscript:
[Detecting RyR clusters with CaCLEAN: influence of spatial distribution and
structural heterogeneity](https://doi.org/10.1101/549683)

Also see the original paper proposing the CaCLEAN algorithm:
[An adaptation of astronomical image processing enables characterization and functional 3D mapping of individual sites of excitation-contraction coupling in rat cardiac muscle](https://doi.org/10.7554/eLife.30425)

## Usage

This set of scripts are designed to perform operations on the results of the 
[CardiacCalcium_FiniteElement](https://github.com/dladd/CardiacCalcium_FiniteElement)
model. Alternatively, results of the finite element models may be obtained
from the author's FigShare archives at [TODO: insert DOI].
The CaCLEAN Evaluator also depends on the
[CaCLEAN](https://github.com/qhtian/CaCLEAN) algorithm. The scripts currently
expect these repositories to be cloned into the same parent directory. Child
dependencies for these projects may also need to be satisfied.

### Confocal Microscopy Simulator

The Python script SimulateMicroscopy.py reads in the finite element nodal values
and interpolates the fluorophore-bound calcium (FCa) field results onto a regular grid
(53nm resolution in all directions)
using discrete natural neighbor (Sibson) interpolation. These results are
convolved with a point spread function to produce optical blurring typical of
confocal microscopy. They are then downsampled to produce 22 images with 215
nm resolution in the image x and y plane and 5 ms temporal resolution. Light noise (SNR=100) is added to the
image data and a background image normalised to the same image resolution and
the FCa initial values. Pixel locations in each image are checked against a
surface model of the cell membrane to construct a boolean mask for the
intracellular space.

The simulated microscopy results and known modeled RyR locations are saved in a Matlab .mat file, to be used by
the CaCLEAN Evaluator.

Python modules used:

* [NumPy](https://www.numpy.org/)
* [SciPy](https://www.scipy.org/)
* [pandas](https://pandas.pydata.org/)
* [naturalneighbor](https://pypi.org/project/naturalneighbor/)
* [trimesh](https://trimsh.org/)
* [skimage](https://scikit-image.org/)
* [math](https://docs.python.org/3/library/math.html)
* [os](https://docs.python.org/3/library/os.html)

### CaCLEAN Evaluator

The TestCaCLEAN_SimulatedMicroscopyResults.m Matlab script operates on the
results of the Confocal Microscopy Simulator using the
[CaCLEAN](https://github.com/qhtian/CaCLEAN) algorithm.

A statistical classification approach was used to assess the performance of RyR cluster detection using CaCLEAN. Modeled cluster centers within the admissible window were considered the actual / ground truth class: TP (ground truth). CaCLEAN detection results were considered the predicted class. Detected RyR cluster sites were defined as determined by the CaCLEAN CRUProps function, which segments cluster regions using Matlab’s built-in watershed algorithm and identifies centroids of segmented regions.

For each modeled cluster location, a TP (detected) classification was assigned if a TP (ground truth) cluster center lied within an available segmented CaCLEAN cluster region (or within a 1 pixel tolerance). When more than one TP (ground truth) fell within a CaCLEAN-detected cluster region, the detected cluster with the nearest centroid to the TP (ground truth) location was marked TP (detected). After classification as TP (detected), the associated CaCLEAN-detected site would be removed from the list of available matches. After iterating through the TP (ground truth) clusters, remaining TP (ground truth) unmatched with detected clusters were classified as false negatives (FN). Remaining CaCLEAN-detected sites unmatched with TP (ground truth) were classified as false positives (FP).

Matlab dependencies used in this project:

* [shadedErrorBar](https://au.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar)
* [colormap](https://au.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-2-0-colormaps-perceptually-uniform-and-beautiful)
* [CaCLEAN](https://github.com/qhtian/CaCLEAN)
* [smoothn](https://au.mathworks.com/matlabcentral/fileexchange/25634-smoothn)(required by CaCLEAN)
* [FMINSEARCHBND](https://au.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon)(required by CaCLEAN)

## Author

* **David Ladd** - [dladd](https://github.com/dladd)

## License

This project is licensed under the Apache 2.0 License - see the LICENSE file for details.

## Acknowledgments

* [Systems Biology Laboratory](https://systemsbiologylaboratory.org.au/) at the University of Melbourne
* [Centre of Excellence in Convergent Bio-Nano Science and Technology](https://www.cbns.org.au/)
* [Cell Structure and Mechanobiology Group](https://cellularsmb.org/) at the University of Melbourne
