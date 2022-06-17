function [ out ] = vascClust (in)
%Requires a structure as input with the following fields
%   fslDir              : top level directory for fsl
%   subject             : string with subject name (e.g., 'FDG001')
%   subOutPath          : subject-level output directory. Created in stareAnchoring.m
%   axialSlices2Clip    : Do not define if you do not want to clip any slices. if desired, the user can elect to collect axial slices from the bottom of the image. Sometimes reconstruction effects are present and this is deriable for clustering
%   midtime             : string of midtimes for all dynamic frames considered of size (n,1) where n is # of midtimes
%   petImagePaths       : cell array of 3D pet image filenames (should be ".hdr" or ".nii" or ".nii.gz")
%   petUnits            : units of pet images. 1 = kBq, 2 = Bq, 3 = mCi

setenv('FSLDIR',in.fslDir);  % tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % tell what the output type should be

fprintf('\n------------------------------------------\nInitiating routine to generate 2 step Vascular Cluster.\n------------------------------------------\n\n')     

out=in;
fprintf('Writing output to: %s.\n',in.subOutPath)
%%%%Create output directeries, convert analyze to nii gz if needed, create 4D
%unzipped nii image
out.vascClust.niiPath=fullfile(in.subOutPath,'anchoring','raw_nii.gz');
if exist(out.vascClust.niiPath,'dir')
    fprintf('Writing to existing nii.gz output directory:\n%s/.\n',out.vascClust.niiPath)
else
    mkdir(out.vascClust.niiPath);
    fprintf('Created and writing to raw moco nii.gz output directory:\n%s/.\n',out.vascClust.niiPath)
end
if length(in.petImagePaths)>1
    if endsWith(in.petImagePaths{1},'.hdr') || endsWith(in.petImagePaths{1},'.nii')
        fprintf('Changing file type from .img or .nii to .nii.gz for all mocoed PET images.\n\n')
        out.vascClust.niiImagePaths={};
        for im=1:length(in.petImagePaths)
            out.vascClust.niiImagePaths{im}=fullfile(out.vascClust.niiPath,[in.subject '.' num2str(im) '.nii.gz']);
            %Convert to nii gz for 3D volumes if not already.
            fsl_cmd=[fullfile(in.fslDir,'bin','fslchfiletype') ' NIFTI_GZ ' in.petImagePaths{im} ' ' out.vascClust.niiImagePaths{im}];
            system(fsl_cmd);
        end
    elseif endsWith(in.petImagePaths{1},'.nii.gz')
        fprintf('File types already .nii.gz for all 3D mocoed PET images. Copying files to appropriate folder.\n\n')
        out.vascClust.niiImagePaths={};
        for im=1:length(in.petImagePaths)
            out.vascClust.niiImagePaths{im}=fullfile(out.vascClust.niiPath,[in.subject '.' num2str(im) '.nii.gz']);
            system(['cp ' in.petImagePaths{im} ' ' out.vascClust.niiImagePaths{im}])
        end
    end
    
    fprintf('Merging all .nii.gz 3D  mocoed files to single 4D image.\n')
    out.vascClust.nii4dPath=fullfile(out.vascClust.niiPath,[in.subject '.nii']);
    %Create 4D nii
    fsl_cmd=[fullfile(in.fslDir,'bin','fslmerge') ' -t ' out.vascClust.nii4dPath ' ' strjoin(out.vascClust.niiImagePaths)];
    system(fsl_cmd);
    system(['gunzip ' out.vascClust.nii4dPath]);    
elseif length(in.petImagePaths)==1
    error('Single PET image provided at in.petImagePaths (either 4D or a single 3D). Please provide a series of 3D images')
end
%%%%

%Load PET image data.
pet4dStruct=load_untouch_nii(out.vascClust.nii4dPath);
%%%%Convert units of PET data to mCi.
if in.petUnits==1
    pet4dData=double(pet4dStruct.img)./37000;
elseif in.petUnits==2
    pet4dData=double(pet4dStruct.img)./37000000;
elseif in.petUnits==3
    pet4dData=double(pet4dStruct.img);
else
    error('in.petUnits set to value other than 1 (KBq), 2 (Bq), or 3 (mCi). Set it correctly or change code in vascClust to accomodate different PET units')
end
%%%%

%%%% Clip bottom axial slices from PET frames if the user elects to do so.
if isfield(in,'axialSlices2Clip') && in.axialSlices2Clip > 0
    fprintf('\nUser selected to clip %d axial slices from bottom of FOV. Presumably because reconstruction can introduce artifacts to the edge slices.\n',in.axialSlices2Clip)
    %Save original PET data structure so that later, we can output clustering
    %images that are the new clipped size, as well as the original size
    %padded with zeros, for ease of viewing the clustering images relative
    %to the original PET data.
    originalPet4dStruct=pet4dStruct;
    %Clip slices & set size of axial dimension to new dimensions
    pet4dData=pet4dData(:,:,in.axialSlices2Clip+1:end,:); 
    pet4dStruct.hdr.dime.dim(4)=size(pet4dData,3);
else
    fprintf('\nUser selected to clip 0 axial slices from bottom of FOV. Presumably because user confirmed that recon edge artifacts were not major factor.\n')
end
%%%%

%Create 2D array of voxel-wise 4D imaging matrix for entry into kmeans
ind=1;
toCluster=zeros(size(pet4dData,1)*size(pet4dData,2)*size(pet4dData,3),size(pet4dData,4));
for d=1:size(pet4dData,3)
    for r=1:size(pet4dData,1)
        for c=1:size(pet4dData,2)
            toCluster(ind,:)=reshape(pet4dData(r,c,d,:),1,size(in.midtime,1));
            ind=ind+1;
        end
    end
end

out.vascClust.figsPath=fullfile(in.subOutPath,'anchoring','figs-masks');
if exist(out.vascClust.figsPath,'dir')
    fprintf('Writing to existing output directory:\n%s/.\n',out.vascClust.figsPath)
else
    mkdir(out.vascClust.figsPath);
    fprintf('Created and writing to output directory:\n%s/.\n',out.vascClust.figsPath)
end

%%%%Run k-means with parallel processing & identify potential vascular clusters in each k-means iteration 
opts=statset('UseParallel',1);
tic
ks=6:4:40; %Running k-means with this range of k clusters
out.vascClust.clustersTested=ks;
kmeansClusIdxs=zeros(size(toCluster,1),length(ks));
out.vascClust.step1MeanTacs=zeros(length(in.midtime),max(ks),length(ks));
out.vascClust.vascClustInds=zeros(1,length(ks)).*NaN;
out.vascClust.vascClustsPeakInd=out.vascClust.vascClustInds;
out.vascClust.vascClustsPeakVal=out.vascClust.vascClustsPeakInd;
for k=1:length(ks)
    fprintf('Working on Step 1 k=%d clusters.\n',ks(k))
    %Run kmeans for given k
    kmeansClusIdxs(:,k)=kmeans(toCluster,ks(k),'emptyaction','drop','MaxIter',1000000,'Replicates',3,'Options',opts);
    %Generate mean TACs from each of the k-means clusters for given k
    for j=1:ks(k)
        out.vascClust.step1MeanTacs(:,j,k)=mean(toCluster(kmeansClusIdxs(:,k)==j,:));
    end
    %Execute the following criteria to identify any potential vascular TAC within all mean TACs for given k (if there are any)
    potentialVascs=out.vascClust.step1MeanTacs(:,:,k); %Grab all of the mean TACs for the given k
    potentialVascs(:,max(potentialVascs)==potentialVascs(end,:))=NaN; %Eliminate TACs where the max is the last time-point (likely irreversible tissue kinetics)
    potentialVascs(:,sum(potentialVascs(2:end,:)<0)>0)=NaN; %Eliminate TACs that have negative values after the 1st time-point (likely noise).
    [tacMaxVal,tacMaxInd]=max(potentialVascs); %tacMaxInd will be 1 for the TACs that are already eliminated and set to NaN.
    tacMaxInd(isnan(potentialVascs(1,:)))=NaN; %In case a TAC that hasn't been eliminated actually has their peak at the 1st time-point, set eliminated TACs back to NaN.
    try
        %In case there are multiple TACs that pass the above tests to be
        %deemed a "vascular TAC", find the vascular TAC with the max early
        %value (most likely indicative of arterial signal, rather than
        %venous or sinus).
        [~,out.vascClust.vascClustInds(k)]=find( max(tacMaxVal(tacMaxInd==min(tacMaxInd))) == tacMaxVal ); %Of the potential vasc TACs that all peak at the earliest time-point, find the TAC with the maximum value.
        %Gather the signal value of the final vascular TAC peak, to be
        %compared with the final vascular TAC peaks from the other
        %k-means iterations.
        [out.vascClust.vascClustsPeakVal(k),out.vascClust.vascClustsPeakInd(k)]=max(out.vascClust.step1MeanTacs(:,out.vascClust.vascClustInds(k),k));
    catch
        fprintf('No vascular TAC identified for k=%d clusters.\n',ks(k))
        out.vascClust.vascClustInds(k)=nan;
        out.vascClust.vascClustsPeakInd(k)=nan;
    end
end
toc
%%%%
if all(isnan(out.vascClust.vascClustsPeakVal))
    error('There was not a single vascular TAC identified in any of the Step 1 cluster options. Quitting. Will need to alter vascular selection criteria or something might be off with the input data.')
end
%%%

%%%Of all of the vascular TACs identified from the different k-means
%%%iterations, find the indice of the most frequent peak & select final
%%%step 1 vascular cluster based on highest peak signal from the frequent
%%%indice.
out.vascClust.frequentPeakInd=mode(out.vascClust.vascClustsPeakInd);
modeInds=out.vascClust.vascClustsPeakInd==out.vascClust.frequentPeakInd;
indOptimalNumClus=find(out.vascClust.vascClustsPeakVal==max(out.vascClust.vascClustsPeakVal(modeInds)));

out.vascClust.step1OptimalNumClus=ks(indOptimalNumClus);
out.vascClust.step1IndOptimalVascClus=out.vascClust.vascClustInds(indOptimalNumClus);
out.vascClust.Step1_idx=kmeansClusIdxs(:,indOptimalNumClus);
%%%

%From the original 2D clustering array, now only include voxels within the
%Step 1 vascular cluster to pass into Step 2 (ignore all non-vascular
%voxel.
step1ToCluster=toCluster.*repmat(out.vascClust.Step1_idx==out.vascClust.step1IndOptimalVascClus,1,size(toCluster,2));
%Save Step 1 Vascular TAC
out.vascClust.step1VascMeanTac=out.vascClust.step1MeanTacs(:,out.vascClust.step1IndOptimalVascClus,indOptimalNumClus);

%%%Output plot of all potential vascular cluster TACs from each k-means
%%%iteration.
plotIn.step1MeanTacs=out.vascClust.step1MeanTacs;
plotIn.midtime=in.midtime;
plotIn.ks=ks;
plotIn.vascClustInds=out.vascClust.vascClustInds;
plotIn.step1VascMeanTac=out.vascClust.step1VascMeanTac;
plotIn.step1OptimalNumClus=out.vascClust.step1OptimalNumClus;
plotIn.figsPath=out.vascClust.figsPath;
plotAnchoringStep1VascTacs( plotIn )
close all
%%%

%%%Properly reshape the k cluster indices from the 2D vector to a 3D image
out.vascClust.step1ClusImg=zeros(size(pet4dData(:,:,:,1)));
ind=1;
for d=1:size(pet4dData,3)
    for r=1:size(pet4dData,1)
        for c=1:size(pet4dData,2)
            out.vascClust.step1ClusImg(r,c,d)=out.vascClust.Step1_idx(ind,1);
            ind=ind+1;
        end
    end
end
%%%

%%%Save k-means cluster masks to nifti file
pet4dStruct.img=out.vascClust.step1ClusImg;
pet4dStruct.hdr.dime.dim(5)=1; %B/c single 3D volume (not 4D time-series as was loaded in).
save_untouch_nii(pet4dStruct,fullfile(out.vascClust.figsPath,['Step1_' num2str(out.vascClust.step1OptimalNumClus) 'clusters_mask.nii']))
%Save mask of just vascular cluster
pet4dStruct.img=double(out.vascClust.step1ClusImg==out.vascClust.step1IndOptimalVascClus);
save_untouch_nii(pet4dStruct,fullfile(out.vascClust.figsPath,['Step1_' num2str(out.vascClust.step1OptimalNumClus) 'clusters_Vasc_only_mask-ind' num2str(out.vascClust.step1IndOptimalVascClus) '.nii']))
%%%
%%%Write out k-means clusters masks in original geometry for easy viewing
%%%if user elected to clip any axial slices
if isfield(in,'axialSlices2Clip')
    paddedImg=zeros(originalPet4dStruct.hdr.dime.dim(2:4));
    paddedImg(:,:,in.axialSlices2Clip+1:end)=out.vascClust.step1ClusImg;
    
    originalPet4dStruct.hdr.dime.dim(5)=1;
    originalPet4dStruct.img=paddedImg;
    save_untouch_nii(originalPet4dStruct,fullfile(out.vascClust.figsPath,['Step1_' num2str(out.vascClust.step1OptimalNumClus) 'clusters_mask_ORIGINAL_geometry-NO_Clipping.nii']));

    paddedImg(:,:,in.axialSlices2Clip+1:end)=out.vascClust.step1ClusImg==out.vascClust.step1IndOptimalVascClus;
    originalPet4dStruct.img=paddedImg;
    save_untouch_nii(originalPet4dStruct,fullfile(out.vascClust.figsPath,['Step1_' num2str(out.vascClust.step1OptimalNumClus) 'clusters_Vasc_only_mask-ind' num2str(out.vascClust.step1IndOptimalVascClus) '-ORIGINAL_geometry-NO_Clipping.nii']));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nStep 2: Run k-means a 2nd time with k=4, but just on the Step 1 extracted vascular cluster.\n\n')     

out.vascClust.step2NumClus=4;
out.vascClust.step2Idx=kmeans(step1ToCluster,out.vascClust.step2NumClus,'emptyaction','drop','MaxIter',1000000,'Options',opts,'Replicates',3);

out.vascClust.step2MeanTacs=zeros(length(in.midtime),out.vascClust.step2NumClus);
for i=1:out.vascClust.step2NumClus
    out.vascClust.step2MeanTacs(:,i)=mean(toCluster(out.vascClust.step2Idx==i,:));
end

%%%As in Step 1, of all of the vascular TACs identified from k=4, find the 
%%%indice of the most frequent peak & select final step 2 vascular cluster 
%%%based on highest peak signal from the frequent indice.
potentialVascs=out.vascClust.step2MeanTacs;
potentialVascs(:,max(potentialVascs)==potentialVascs(end,:))=NaN; %Eliminate TACs where the maximum value is at the end-point
potentialVascs(:,sum(potentialVascs(2:end,:)<0)>0)=NaN; %Eliminate TACs where anything but the 1st point is less than zero
%Find maximum value of each TAC
[tacMaxVal,tacMaxInd]=max(potentialVascs); %tacMaxInd will be 1 for the TACs that are already eliminated and set to NaN.
tacMaxInd(isnan(potentialVascs(1,:)))=NaN; %In case a TAC that hasn't been eliminated actually has their peak at the 1st time-point, set NaN'ed TACs back to NaN.
try
    [~,out.vascClust.step2IndOptimalVascClus]=max(tacMaxVal); %find the TAC with the maximum value

    out.vascClust.step2VascMeanTac=out.vascClust.step2MeanTacs(:,out.vascClust.step2IndOptimalVascClus);
catch
    warning('Step 2 k-means did not yield IDIF meeting appropriate criteria.')
    save(fullfile(in.subOutPath,'anchoring',[in.subject '_anchoring_out.vascClust.mat']),'out')
    keyboard
end
%%%

%%%Properly reshape the k cluster indices from the 2D vector to a 3D image
out.vascClust.step2ClusImg=zeros(size(pet4dData(:,:,:,1)));
ind=1;
for d=1:size(pet4dData,3)
    for r=1:size(pet4dData,1)
        for c=1:size(pet4dData,2)
            out.vascClust.step2ClusImg(r,c,d)=out.vascClust.step2Idx(ind,1);
            ind=ind+1;
        end
    end
end
%%%

%%%Plot final 2-step k-means vascular TACs
plotIn={};
plotIn.step1VascMeanTac=out.vascClust.step1VascMeanTac;
plotIn.step2VascMeanTac=out.vascClust.step2VascMeanTac;
plotIn.midtime=in.midtime;
plotIn.figsPath=out.vascClust.figsPath;
plotIn.subject=in.subject;
plotAnchoringStep2Step1VascTacs( plotIn )
close all
%%%

%%%Save k-means cluster masks to nifti file
pet4dStruct.img=out.vascClust.step2ClusImg;
save_untouch_nii(pet4dStruct,fullfile(out.vascClust.figsPath,['Step2_' num2str(out.vascClust.step2NumClus) 'clusters_mask.nii']))
%Save mask of just vascular cluster
pet4dStruct.img=double(out.vascClust.step2ClusImg==out.vascClust.step2IndOptimalVascClus);
out.vascClust.step2VascMaskPath=fullfile(out.vascClust.figsPath,['Step2_' num2str(out.vascClust.step2NumClus) 'clusters_Vasc_only_mask-ind' num2str(out.vascClust.step2IndOptimalVascClus) '.nii']);
save_untouch_nii(pet4dStruct,out.vascClust.step2VascMaskPath)
%%%
%%%Write out k-means clusters masks in original geometry for easy viewing
%%%if user elected to clip any axial slices
if isfield(in,'axialSlices2Clip')
    paddedImg=zeros(originalPet4dStruct.hdr.dime.dim(2:4));
    paddedImg(:,:,in.axialSlices2Clip+1:end)=out.vascClust.step2ClusImg;
    
    originalPet4dStruct.hdr.dime.dim(5)=1;
    originalPet4dStruct.img=paddedImg;
    save_untouch_nii(originalPet4dStruct,fullfile(out.vascClust.figsPath,['Step2_' num2str(out.vascClust.step2NumClus) 'clusters_mask_ORIGINAL_geometry-NO_Clipping.nii']));

    paddedImg(:,:,in.axialSlices2Clip+1:end)=out.vascClust.step2ClusImg==out.vascClust.step2IndOptimalVascClus;
    originalPet4dStruct.img=paddedImg;
    out.vascClust.step2VascMaskOrigGeomPath=fullfile(out.vascClust.figsPath,['Step2_' num2str(out.vascClust.step2NumClus) 'clusters_Vasc_only_mask-ind' num2str(out.vascClust.step2IndOptimalVascClus) '-ORIGINAL_geometry-NO_Clipping.nii']);
    save_untouch_nii(originalPet4dStruct,out.vascClust.step2VascMaskOrigGeomPath);
end
delete(gcp('nocreate'))
save(fullfile(in.subOutPath,[in.subject '_STARE_out.mat']),'out')

fprintf('\n------------------------------------------\nDONE. Vascular cluster extraction complete.\n------------------------------------------\n\n')     
