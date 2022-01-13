function [ stareOut ] = stare( stareIn )
%stare Runs all necessary functions for Source-to-Target Rotating
%Estimation.
%For information on the algorithm, see Bartlett et al. 2022 in NeuroImage.
%   Requires a structure as input with the following fields:
%   outPath             : String that is top level directory where you want subject-specific STARE output folder to be generated.
%   subject             : String with subject name (e.g., 'FDG001')
%   axialSlices2Clip    : Do not define if you do not want to clip any slices. if desired, the user can elect to clip axial slices from the bottom of the PET images. Sometimes reconstruction effects are present in the inferior-most slices and this is undesirable for whole field of view (FOV) clustering. 6 slices were clipped in the filtered back projection (FBP) images in Bartlett et al. 2022
%   midtime             : String of midtimes for all dynamic frames considered of size (m,1) where m is # of midtimes
%   fslDir              : Top level directory for fsl on your machine (e.g., '/usr/local/fsl/')
%   petImagePaths       : Cell array of strings taht are the 3D pet image filenames (should be ".hdr" or ".nii" or ".nii.gz")
%   petUnits            : Units of pet images. 1 = kBq, 2 = Bq, 3 = mCi
%   pvcMethod           : Must be 'STC'. The PETPVC toolbox has many PVC options. STARE has only been validated with 'STC' - Single Target Correction.
%   fwmh                : The full width at half maximum (fwhm) of the point spread function of the PET images. This should account for scanner resolution and reconstruction approach. For details see https://github.com/UCL/PETPVC and Sari et al. 2017 "Estimation of an image derived input function with MR-defined carotid arteries in FDG-PET human studies using a novel partial volume correction method".
%   regions             : A cell array of strings with the brain region names that will be quantified in STARE.
%   tacs                : Time activity curves (TACs) for brain regions in "regions" (in the same order) as m x r array where m is length(midtime) and r is length(regions)
%   tracer              : Radiotracer used. Only validated for FDG currently.
%   vascCorrPerc        : Percentage from 0 to 100 (0 and 5 used in validaions in Bartlett et al. 2022) for contribution to voxel-wise PET signal from vasculature signal (i.e., if =5, 5% of signal is corrupted from vasclature and vascular correction of TACs will be run with Vb = 5.

fprintf('Running STARE for subject: %s.\n',stareIn.subject)  
%From the top level user-defined output directory, create a
%subject-specific output directory with time info
stareIn.subOutPath=fullfile(stareIn.outPath,[stareIn.subject '_' datestr(now,'yyyymmdd_HHMM')]);
if exist(stareIn.subOutPath,'dir')
    fprintf('Writing to existing output directory:\n%s/.\n',stareIn.subOutPath)
else
    mkdir(stareIn.subOutPath);
    fprintf('Created and writing to output directory:\n%s/.\n',stareIn.subOutPath)
end
%%%%%%%%%%%%%%%%%%%%%%%%% Now run through all STARE modules, passing the
%%%%%%%%%%%%%%%%%%%%%%%%% output of one as the input to the next.

%Run vascular clustering
[ vcOut ] = vascClust(stareIn);

%Run partial volume correction of vascular cluster.
[ petPvcOut ] = petPvc(vcOut);

%Use output of petPvc to vascularly correct TACs (by extracting the mean
%signal in the PVC'ed vasculature cluster, fitting it, and then applying
%vascular correction
[ fitVcOut ] = fitVascMeanTac(petPvcOut);
[ vascCorrOut ] = tacVascCorr(fitVcOut);
%Bootstrap signal in PVC'ed vasculature to generate many potential input
%functions to compartmental modeling, that will then allow us to derive
%ranges of rate constants and the penalty term to be used in the STARE cost
%function
[ bootOut ] = bootAnchor(vascCorrOut);
%Now that we have the subject-specific anchoring, minimize the cost
%function
[ stareOut ] = minimizeCostFunc(bootOut);














