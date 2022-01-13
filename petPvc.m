function [ out ] = petPvc(in)
%petPvc Formats files and runs PETPVC toolbox on the final vasculature
%cluster identified in vascClustering.
%   To be run after vascClustering. Can pass output structure of
%   vascClustering as input to petPvc and add the following 3 additional
%   fields, or assign each field required individually as shown below.
%
%   vascClustering output structure with these 3 fields added.
%   axialSlices2Clip        : Do not define if you do not want to clip any slices. if desired, the user can elect to collect axial slices from the bottom of the image. Sometimes reconstruction effects are present and this is deriable for clustering
%   pvcMethod               : Must be 'STC'. The PETPVC toolbox has many PVC options. STARE has only been validated with 'STC' - Single Target Correction.
%   fwmh                    : float for the full width at half maximum (fwhm) of the point spread function of the PET images. This should account for scanner resolution and reconstruction approach. For details see https://github.com/UCL/PETPVC and Sari et al. 2017 "Estimation of an image derived input function with MR-defined carotid arteries in FDG-PET human studies using a novel partial volume correction method".
% OR
%   The above 3 fields plus each of the following fields:
%   fslDir                  : top level directory for fsl
%   subOutPath              : subject-level output directory. Created in stareAnchoring.m
%   If axialSlices2Clip is empty (don't want to clip any slices)
%        vascClust.step2VascMaskPath              : 3D mask with zeros in background and ones in the locations of the vasculature cluster (in original geometry/dimensions of the clipped PET data)
%   If axialSlices2Clip is set to value (want to clip any slices)
%        vascClust.step2VascMaskOrigGeomPath      : 3D mask with zeros in background and ones in the locations of the vasculature cluster (in original geometry/dimensions of the PET data)
%   vascClust.niiImagePaths           : a cell array of full file paths to each of the 3D nii PET images (entire time-series). petpvc requires 3D files, rather than a 4D nii of the entire scan
%   petUnits                : units of pet images. 1 = kBq, 2 = Bq, 3 = mCi
%   vascClust.step1VascMeanTac        : time activity curve (TAC) from the final step 1 vascular cluster
%   VascClust.step2VascMeanTac        : time activity curve (TAC) from the final step 2 vascular cluster
%   midtime                 : vector of midtimes for all dynamic frames considered of size (n,1) where n is # of midtimes
%   VascClust.figsPath                : ouput path
%   subject                 : subject ID

fprintf('\n------------------------------------------\nInitiating routine to Partial Volume Correct (PVC) 2-step vascular cluster previously generated.\n------------------------------------------\n\n')     

setenv('FSLDIR',in.fslDir);  % tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % tell what the output type should be

%Check to ensure pvcMethod and fwhm are set.
if in.pvcMethod~='STC'
    error('Please set in.pvcMethod to STC. PETPVC toolbox has other options for pvc, but these are not currently added in STARE.')
end
if ~isfield(in,'fwhm')
    error('Please set a float value for: fwhm.')
end

out=in;
%Create PVC output directory.
out.petPvc.pvcOutPath=fullfile(in.subOutPath,'anchoring','PVC');
if exist(out.petPvc.pvcOutPath,'dir')
    fprintf('Writing to existing output directory:\n%s/.\n\n',out.petPvc.pvcOutPath)
else
    mkdir(out.petPvc.pvcOutPath);
    fprintf('Created and now writing to output directory:\n%s/.\n\n',out.petPvc.pvcOutPath)
end
%Identify final vascular cluster mask
if isfield(in,'axialSlices2Clip')
    vascMaskPath=in.vascClust.step2VascMaskOrigGeomPath;
else
    vascMaskPath=in.vascClust.step2VascMaskPath;
end
        
fprintf('Using following file as vascular cluster mask:\n%s\n',vascMaskPath)
%%%Run PVC
fprintf('Running petpvc algorithm on vascular cluster masked PET time-series.\n')
out.petPvc.niiPvcOut={};
for im=1:length(in.vascClust.niiImagePaths)
    [~,niiName,niiExt]=fileparts(in.vascClust.niiImagePaths{im});
    out.petPvc.niiPvcOut{im}=fullfile(out.petPvc.pvcOutPath,['PVC_' niiName niiExt]);
    petPvcCmd=['/Users/Betsy/Desktop/PETPVC-build/src/petpvc -i ' in.vascClust.niiImagePaths{im} ' -o ' out.petPvc.niiPvcOut{im} ' -m ' vascMaskPath ' -p ' in.pvcMethod ' -x ' num2str(in.fwhm) ' -y ' num2str(in.fwhm) ' -z ' num2str(in.fwhm)];
    if im==1
        fprintf('Sample petpvc command: %s\n',petPvcCmd)
    end
    system(petPvcCmd);
end
%%%
%Write 4D version of PVCed PET images
out.petPvc.pvc4dImage=strrep(out.petPvc.niiPvcOut{im},['.' num2str(im) '.'],'.');
fsl_cmd=[fullfile(in.fslDir,'bin','fslmerge') ' -t ' out.petPvc.pvc4dImage ' ' strjoin(out.petPvc.niiPvcOut)];
system(fsl_cmd);
system(['gunzip ' out.petPvc.pvc4dImage]);

pvc4dStruct=load_untouch_nii(strrep(out.petPvc.pvc4dImage,'.gz',''));
%%%%Convert units of PET data to mCi.
if in.petUnits==1
    pvc4dData=double(pvc4dStruct.img)./37000;
elseif in.petUnits==2
    pvc4dData=double(pvc4dStruct.img)./37000000;
elseif in.petUnits==3
    pvc4dData=double(pvc4dStruct.img);
else
    error('in.petUnits set to value other than 1 (KBq), 2 (Bq), or 3 (mCi). Set it correctly or change code in vascClustering to accomodate different PET units')
end
%Load in mask for final vascular cluster
vascMaskStruct=load_untouch_nii(vascMaskPath);
vascMaskData=double(vascMaskStruct.img);

reshapedPvc4dData=reshape(pvc4dData,size(pvc4dData,1)*size(pvc4dData,2)*size(pvc4dData,3),size(pvc4dData,4));
%Get TACs of PVCed PET data in the final vascular cluster
out.petPvc.pvcVascMeanTac=mean(reshapedPvc4dData(vascMaskData==1,:))';
out.petPvc.pvcVascSdTac=std(reshapedPvc4dData(vascMaskData==1,:))';

%%%Plot final PVCed vascular TACs
plotIn={};
plotIn.step1VascMeanTac=in.vascClust.step1VascMeanTac;
plotIn.step2VascMeanTac=in.vascClust.step2VascMeanTac;
plotIn.pvcVascMeanTac=out.petPvc.pvcVascMeanTac;
plotIn.midtime=in.midtime;
plotIn.figsPath=in.vascClust.figsPath;
plotIn.subject=in.subject;
plotAnchoringPvcVascTacs( plotIn )
close all
%%%

save(fullfile(in.subOutPath,[in.subject '_STARE_out.mat']),'out')

fprintf('\n------------------------------------------\nDONE. 2 step vascular cluster now PVCed and TACs extracted.\n------------------------------------------\n\n')     
