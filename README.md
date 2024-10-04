# Luminance-Contrast-Aware Foveated Rendering

This repository contains the Matlab implementation of the sigma predictor from:

Tursun, O. T., Arabadzhiyska-Koleva, E., Wernikowski, M., Mantiuk, 
R., Seidel, H. P., Myszkowski, K., & Didyk, P. (2019). 
Luminance-contrast-aware foveated rendering. ACM Transactions on 
Graphics (TOG), 38(4), 98.

In order to run the predictor please download matlabPyrTools by Simoncelli
from 

    https://github.com/LabForComputationalVision/matlabPyrTools

and place it in the folder: `./matlabpyrtools`

Then compile the MEX-files by running the following script:

    ./matlabpyrtools/MEX/compilePyrTools.m

and place the created files (\*.mexw64, \*.mexa64 or \*.mexmaci64 depending
on the platform) in Matlab path (see path, addpath, genpath and pathtool
functions of Matlab).

# Files
`main.m` : Provides the implementation of sigma predictor for an image patch.
Feel free to use "help main" for info.

`run_on_image.m`: Runs the predictor on patches of an image and returns the
predicted sigma map.

`sample_run.m`: Loads a sample image from the paper and runs the predictor.
Try running this first.

`get_params.m`: Loads optimum predictor parameters learned using Simulated
Annealing.

`display_params.m`: Returns the parameters such as resolution, physical size
and the observation distance for the particular display used in our
experiments. If the physical parameters of your display are significantly
different from ours, you will need to define your display in this file and
modify `run_on_image.m:11-31` where those parameters are used.

`disp1_luminance.mat` and `disp2_luminance.mat`: Provides the calibration
data for converting RGB to luminance for the displays that we used in our
experiments. The predictor expects luminance (cd/m^2) as the input. If 
the calibration data does not exist, inverse gamma transformation may be 
used to approximately compute linear values from sRGB (gamma = 2.2).

# OpenGL and Unity Applications

They are available for download at:

    https://www.pdf.inf.usi.ch/projects/AdaptiveFoveation/AdaptiveFoveation.zip

# Links

Our project website:

    https://www.pdf.inf.usi.ch/projects/AdaptiveFoveation

