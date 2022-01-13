STARE - Source-to-Target Rotating Estimation

https://www.sciencedirect.com/science/article/pii/S1053811922000313

STARE has been developed and validated in Matlab R2016b.

Software dependencies:
FSL 6.0.01
NifTI 2014 - included in subfolder to STARE code (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image?s_tid=prof_contriblnk)

Add STARE code folder and subolders to Matlab path first.

Start with wrapperStare.m.

STARE requires the PET midtimes, 3D files of the dynamic PET data, and time activity curves for the regions for which quantification is sought. Examples of pulling this data into matlab and passing it into STARE in the correct format is provided in wrapperStare.m

STARE also requires a set of user-defined inputs. These are defined and a sample set of user inputs are also provided in wrapperStare.m. Start here. wrapperStare.m calls stare.m, which is the main function and contains a set of modularized functions to run STARE from start to finish.

