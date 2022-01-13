%%% A wrapper for STARE
%For information on the algorithm, see Bartlett et al. 2022 in NeuroImage.

%%stare.m will be the function run (last lines of this script) and its
%%inputs, as defined at the top of the stare.m function are also copied
%%here for convenience and defined with examples in this script.

%To run through all components of STARE, the following inputs are required.
%   Requires a structure as input with the following fields:
%   outPath             : String that is top level directory where you want subject-specific STARE output folder to be generated.
%   subject             : String with subject name (e.g., 'FDG001')
%   axialSlices2Clip    : Do not define if you do not want to clip any slices. if desired, the user can elect to clip axial slices from the bottom of the PET images. Sometimes reconstruction effects are present in the inferior-most slices and this is undesirable for whole field of view (FOV) clustering. 6 slices were clipped in the filtered back projection (FBP) images in Bartlett et al. 2022
%   midtime             : String of midtimes for all dynamic frames considered of size (m,1) where m is # of midtimes
%   fslDir              : Top level directory for fsl on your machine (e.g., '/usr/local/fsl/')
%   petImagePaths       : Cell array of strings that are the motion-corrected 3D pet image filenames (should be ".hdr" or ".nii" or ".nii.gz")
%   petUnits            : Units of pet images. 1 = kBq, 2 = Bq, 3 = mCi
%   pvcMethod           : Must be 'STC'. The PETPVC toolbox has many PVC options. STARE has only been validated with 'STC' - Single Target Correction.
%   fwmh                : The full width at half maximum (fwhm) in millimeters of the point spread function of the PET images. This should account for scanner resolution and reconstruction approach. For details see https://github.com/UCL/PETPVC and Sari et al. 2017 "Estimation of an image derived input function with MR-defined carotid arteries in FDG-PET human studies using a novel partial volume correction method".
%   regions             : A cell array of strings with the brain region names that will be quantified in STARE.
%   tacs                : Time activity curves (TACs) for brain regions in "regions" (in the same order) as m x r array where m is length(midtime) and r is length(regions)
%   tracer              : Radiotracer used. Only validated for FDG currently.
%   vascCorrPerc        : Percentage from 0 to 100 (0 and 5 used in validaions in Bartlett et al. 2022) for contribution to voxel-wise PET signal from vasculature signal (i.e., if =5, 5% of signal is corrupted from vasclature and vascular correction of TACs will be run with Vb = 5.

%Below entries for each of the above user defined inputs are just examples
%and need to be replaced with real user inputs to be able to run algorithm.
stareIn={};
stareIn.outPath='/tmp/stareOut/';
stareIn.subject='FDG001';
stareIn.axialSlices2Clip=6;
stareIn.midtime=load(fullfile('/tmp/stareIn/',stareIn.subject,'midtimes.txt'));
stareIn.fslDir='/usr/local/fsl';

imagePath=fullfile('/tmp/stareIn/',stareIn.subject,'/mocoPetData/');
imgDir=dir(fullfile(imagePath,'*.nii'));
stareIn.petImagePaths=fullfile(imagePath,{imgDir.name})';

stareIn.petUnits=1; %PET images in Bartlett et al. 2022 were kBq

stareIn.pvcMethod = 'STC'; %The PETPVC toolbox has many PVC options. STARE has only been validated with 'STC' - Single Target Correction
stareIn.fwhm=5.9; %PET images used in Bartlett et al. 2022 required a fwhm of 5.9mm

stareIn.regions={'cerfullcs_c','cin','hip','par','pfc','pip'}; %These were the regions validated in Bartlett et al. 2022. Switch to your own atlas region names

%Extract TAC data from all regions
stareIn.tacs=importdata(fullfile('/tmp/stareIn/',stareIn.subject,'tacs.txt'));

stareIn.tracer='FDG';
stareIn.vascCorrPerc=5;

%Ready to run!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[stareOut] = stare(stareIn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%