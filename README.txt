CALCIUM
DOI: 10.5281/zenodo.6390856

OVERVIEW:

This code analyses cell movement and intracellular calcium release. It’s designed for fluorescence microscopy experiments where cells are labelled with a calcium indicator and dropped onto a surface (e.g. antibody-labelled glass or a supported lipid bilayer). For cells which are stationary throughout (and in place from the start) use calciumStationaryCells.m instead.

INSTRUCTIONS:

Change the input parameters at the top of the main file (calcium.m). Make sure the input directory contains a single TIFF file.

DESCRIPTION:

The code reads in a single TIFF file representing a 3D time series. After background subtraction, flat-field correction and Gaussian smoothing, local maxima in each frame are identified. These are recorded (for a subset of frames) in Peaks.tif.

Maxima from each frame are combined into tracks using a nearest-neighbour approach. Tracks shorter than a user-defined length are discarded. Remaining tracks are recorded in Cells.tif and Labelled_cells.tif.

The distance travelled by each cell over time is calculated. A user-defined speed threshold is used to identify the time at which cells become attached to the substrate. To be classified as attached a cell must remain beneath the threshold speed until the end of the track. Graphs of distance and speed with respect to time are saved in the Distance folder.

The mean intensity of a cell-sized disk around each peak is calculated, forming intensity traces for each cell. Noise is reduced by Gaussian smoothing. Sharp increases in intensity are identified from the first derivative and a user-defined threshold determines which of these are classified as calcium spikes. If a minimum spike intensity and/or minimum spike duration have been specified at the start, spikes not meeting these criteria are discarded.

Cells are classified as “coming in hot” when no spikes are indentified but the intensity at the start of the track exceeds a user-defined threshold. These cells are excluded when calculating the % of cells releasing calcium (“triggering”). To include cells with high intensities at the start among those not releasing calcium, set this threshold to Inf.

The time since the cell was first seen is recorded for each spike, as well as the time since attachment (this will be negative if the cell is still moving). A spike starts at a peak in the first derivative and ends when the intensity drops below its level at the start of the spike. A second spike can only begin once the first has ended. Graphs of intensity and intensity gradient with respect to time are saved in the Calcium_release, No_calcium_release  and Came_in_hot folders.

The duration, maximum intensity and integrated intensity of each spike are recorded. The areas corresponding to integrated intensity are shaded on the intensity/time graphs in grey. The time of attachment is marked with a black diamond. Data relating to individual spikes are recorded in Spikes.xlsx and data relating to individual cells in Results.xlsx. The % of attached/triggered cells is plotted with respect to time in Landing.png/Triggered.png.

AUTHOR:

Jane Humphrey (janehumphrey@outlook.com)

LICENSE:

MIT (see LICENSE.txt file for details)

DEPENDENCIES:

MATLAB v9.5
Signal Processing Toolbox v8.1
Image Processing Toolbox v10.3
Statistics and Machine Learning Toolbox v11.4

The code has been tested on Windows 8, Windows 10 and MacOS operating systems.

INSTALLATION:

MATLAB can be obtained from https://www.mathworks.com/products/matlab.html. No further installation is required.

TEST DATA:

An example TIFF stack can be found at https://figshare.com/articles/Test_data/12298802.

ACKNOWLEDGEMENTS:

The functions pkfnd.m and track.m were written by Daniel Blair and Eric Dufresne. More information can be found at http://site.physics.georgetown.edu/matlab.

The export_fig toolbox was written by Yair Altman (https://github.com/altmany/export_fig).

Thanks to Dr Aleks Ponjavic for inspiration and Prof. Sir David Klenerman for support.
