# dm3ReaderAddon

Standard Operating Procedure: dm3ReaderAddon.m
Author: Phong Hien Nguyen
Date: 02 August 2017
 
Purpose

The intended purpose of this MATLAB program is to extract element-concentration information from Electron Energy Loss Spectrum (EELS) images. 

 
Quick-Start Guide

	Open dm3ReaderAddon.m
	Press F5 to run the program
	Specify the settings to be used (see the Command Window)
	Directory navigation user-interface (UI) will open. Navigate to the high-loss EELS image and select it.
	If a low-loss image was specified available, select it in the next directory navigation UI.
	Enter user-defined parameters (element-specific background estimation and integration parameters).
	Output datasets are typically named after the user-specified element along with “_EELS.csv” and are in the directory that dm3ReaderAddon.m is saved in, in a folder called, ‘dm3ReaderAddon_output’
	The “BRPixels_” prefix denotes that the *.csv file contains background-subtracted EELS data per pixel. The columns, from left-to-right, represent pixels from right-to-left, bottom-to-top. The first column on the left contains the eV-axis that corresponds to every intensity-vector on the right.
	The “DataSquare_” prefix denotes that the *.csv file contains background-subtracted EELS data, summed across pixel-rows (to give top-to-bottom information) or summed across pixel-columns (to give left-to-right information). The first column on the left contains the eV-axis that corresponds to every intensity-vector on the right. 
	The “IntVector_” or “IntVectorNorm” prefix denotes that the *.csv file contains the intensity of the user-specified element. The first column on the left contains the distance-axis that corresponds to every intensity-vector on the right. 
	The “IntSquare_” prefix denotes that the *.csv file contains the intensity of the user-specified element with each intensity maintaining the same coordinate as the pixel in the EELS image if the direction setting is 0. 
 
Requirements and Settings

This program requires the following resources to be present in the same directory as the program file to run:
	dm3ReaderModified.m
	a modified version of dm3Reader.m which outputs a text file containing metadata tied to the *.dm3 file input
	EELSData.m
	a class definition file that defines the structure of objects used in dm3ReaderAddon.m

Alongside those files listed above, user defined inputs are as follows:

	direction
	if top-to-bottom (option 1) is specified, the program will sum background-subtracted, element-integrated values across pixel rows of the image, outputting vertical-distance versus normalized element-concentration data
	if left-to-right (option 0) is specified, the program will sum background-subtracted, element-integrated values across pixel columns of the image, outputting horizontal-distance versus normalized element-concentration data
	LLExist
	if the user specifies that a low-loss spectrum image is available (option 1), the user will need to provide two *.dm3 files: the first is the high-loss EELS image and the second is the low-loss EELS image. The low-loss EELS image is used to perform a zero-loss peak shifting on the corresponding high-loss EELS image pixels. Otherwise (option 0), no zero-loss peak shift is performed.
	BRExist
	if the user specifies that the background should be removed (option 1 or 2), a range for of values to be used for artificial fitting of the background must be provided later. Otherwise (option 0), the background will not be removed.
	Option 1: Power law model: f(x)=bx^m 
	Option 2: Exponential (1st-order log polynomial) model: f(x)=be^mx
	userDefAxis
	if the user specifies that the spatial axis should be manually created (option 1) (such as when spatial drift correction was not used when the EELS image was acquired, or if the spatial drift correction is unreliable), then the user must later specify the physical height of the EELS image (most imaging software will provide tools that will allow a user to measure a scale bar and obtain a distance/pixel ratio. This information, along with image dimensions can be used to figure out what the physical height portrayed in the EELS image is, provided an ADF image with a scale bar was taken concurrently). Output distance-vectors will have the same units as provided for this query; however, the input should not contain any non-numerical inputs.
	Otherwise (option 0), this program will use metadata from the tag file created by dm3ReaderModified.m to construct spatial axis data. This data comes directly from the *.dm3 file and thus, is as accurate as the microscope output (i.e., dependent on accuracy of spatial drift correction).
	norm2Material
	if the user specifies that element counts should be normalized to those within the range of the material (option 1), then all element count data points will be subtracted by the minimum value, then divided by the maximum intensity value within the user-specified range. Though a plot will be generated in an attempt to show the user a preliminary graph of the element distribution as a function of distance, it is recommended that the user know in advance which pixel (top-to-bottom or left-to-right appropriately) the material of interest begins and ends in the EELS image.

 
Main Function

In sequential order, the main function will:
	Ask for the settings the user wishes to use.
	If previously saved settings are available, the user will be prompted to select the settings file. Otherwise, the user will have to enter the user-defined settings and address an option to save these settings at the end.
	Create empty object structures defined by the EELSData class definition
	Request the high-loss EELS image *.dm3 file
	Store the name of the *.dm3 file for later use
	Store the name of the text file containing the corresponding metadata
	Read useful metadata from the tags file:
	HLinput.scale(1) contains the physical size represented by a pixel in the vertical direction (or horizontal direction, depending on whether the user specified top-to-bottom or left-to-right information as desired)
	HLinput.scale(2) contains the physical size represented by a pixel in the direction not addressed by HLinput.scale(1)
	HLinput.scale(3) contains the eV step size used when taking the EELS data for any pixel in the EELS image
	HLinput.origin contains information about the starting point for the EELS data. The true starting point of the EELS data is the product of the (origin*-1) and HLinput.scale(3)
	Request the low-loss EELS image *.dm3 file and read the useful metadata from the tags file
	LLinput.scale(3) contains the eV step size used when taking the EELS data for any pixel in the EELS range
	LLinput.origin contains information about the starting point for the EELS data. The true starting point of the EELs data is the product of the (origin*-1) and LLinput.scale(3)
	Request the name of the first element to be analyzed within the EELS image
	The name of the output files are generated by taking this string input and concatenating the input with “_EELS.csv”. The output CSV file is located in the directory containing dm3ReaderAddon.m and other resources, in a folder named, “dm3ReaderAddon_output”.
	A copy of HLinput is made, called HLoutput
	HLoutput is meant to be used by the program as an adjustable set of data. HLinput is meant to stay pristine, should the raw data be needed later in the program.
	If the user specified that the background should be removed, LLEst and ULEst are requested.
	LLEst and ULEst specify the lower and upper limits of the range of eV values to use for background-estimation. This should be done in advance by fitting background to immediately before the peak of interest (the peak corresponding to the element specified earlier) using Gatan Digital Micrograph.
	eVExtract is called upon. eVExtract takes origin and scale data from the input structure and creates an initial eV vector for the 2048 intensity values in the EELS data for each pixel. At this point, all pixels share the same eV vector.
	backgroundRemove is called upon. For each pixel, artificial eV vectors and intensity vectors are created, using only values from LLEst to ULEst. A power law or decaying exponential model is fit to the artificial vectors. An intensity value is approximated for every point in the full eV vector. This background-approximating vector is then subtracted from the pixel’s actual spectrum and output into HLoutput.Data, resulting in background-removed EELS data. The background-fitting function can be improved upon by modifying it within the backgroundRemove function near the end of the program.
	If a low-loss EELS image was provided, the program will find the zero-loss peaks by querying the data for the highest value and store the eV position of the highest value (if the highest value is shared among multiple eV values, the average will be used). This zero-loss peak-finding function can be improved by taking the average of all eV points within a fit Gaussian’s full width at half maximum.
	After finding the least shift (the peak farthest to the left) and the greatest shift (the peak farthest to the right), the eV axis is extended by defining the first point to be the difference between the first value of the eV axis and the greatest shift. The last value is defined as the difference between the last value and the least shift. The step size between this first and last value is the difference between elements 2 and 1 within the original eV vector.
	To make the EELS data line up with this new unified eV vector, the actual EELS data are padded on the left and right of each intensity vector. The amount of 0s to be padded on the left are defined as the integer step size difference between the shift corresponding to that specific pixel (the specific shift is in output.shiftMatrix, at the same coordinates as the pixel in the EELS image) and the minimum shift of all shifts (this is the distance from the starting eV value in the unified eV axis that our pixel-specific spectrum actually starts). The number of 0s to be padded on the right are defined as the integer step size difference between the maximum shift of all shifts and the shift corresponding to a specific pixel. 
	The integration range is shared amount all spectra in the EELS image. It is defined from user-specified lower limit to the user-specified lower limit + 20 eV. This can be edited directly in this portion of the code. Specifically, a user-defined upper integration limit can be implemented by replacing:
ULInt=LLInt+20;
with 
input('input('UL for integration (multiple of %.3f): \n',HLinput.scale(1));
	Each pixel’s EELS data is integrated using the specified upper and lower limits of integration using a trapezoidal method (built into MATLAB). The specific integration value for an EELS pixel is saved in a matrix with the same dimensions as the image itself, and is called output.intSquare.
	Data is summed across rows to produce information as a function of only row position (rows and columns are switched appropriately depending on if the user specified information pertaining to top-to-bottom or left-to-right). This information is normalized to 1 and the output is placed in output.intVector.
	If the user specifies that the spatial data should be created manually, the user will be prompted for the physical height represented in the EELS image. This information, along with the indices of the image are used to create a spatial vector to put alongside the integrated values (the normalized concentration of a specific element. Otherwise, metadata extracted from the tag file outputted by dm3ReaderModified.m will be used). 
	The purpose of normalizing data to within the material region of the image is to maximize the contrast of information within the material. If the user enables this feature, the data will be normalized such that only the data within the material-region (specified by the user) is considered when normalizing all data to 1. For convenient comparison to EDX images, these coordinates (their real-space equivalents) are stored in output.materialCoordinates with the first and second elements corresponding to the material-start and material-end regions respectively.
	Lastly, data is extracted and quick plots can be enabled. Within this last region of the main function, different structures can be written into *.csv files. The current scheme is to append the name of the structure to the front of the name of the element and “_EELS.csv” to use as the name for the output *.csv file. All structures in the “output” structure can be extracted to *.csv files in this fashion. Furthermore, these writing functions can be turned off with simple if loops, demonstrated with examples in this section of the code. Plots for initial verification of data handling by the code can also be turned on this way. Examples of this are in present within  	 this section as well. 
