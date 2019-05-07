% extract features


%% perturbed TCGA-TCIA data
workingDir = '/Users/messi/Documents/private_research/Stanford_BMI/HeadNeck_Cancer/data/original/hnscc_CT/hnscc_out_imgs/perturbed/';
cd (workingDir);
perturb_paths = dir(workingDir);
% delete '.', '..', '.DS_Store', and '_DS_Store'
toRemoveIdx = [];
for i = 1:length(perturb_paths)
    if startsWith(perturb_paths(i).name,'.') || startsWith(perturb_paths(i).name,'..') || strcmp(perturb_paths(i).name,'.DS_Store') || strcmp(perturb_paths(i).name,'_DS_Store') || strcmp(perturb_paths(i).name,'__.DS_Store') || strcmp(perturb_paths(i).name,'._.DS_Store') || contains(perturb_paths(i).name,'.csv') || contains(perturb_paths(i).name,'.xml') || contains(perturb_paths(i).name,'.mat') % '._.DS_Store' required for Linux.
        toRemoveIdx = [toRemoveIdx,i];
    end;
end;
perturb_paths(toRemoveIdx,:) = [];

outputDir = workingDir;
renew_file = ''; % file to be appended with feature rows for new patients
for i = 1:length(perturb_paths)
    relative_path_to_ROI = perturb_paths(i).name;
    prefix_of_output_csv = perturb_paths(i).name;
    [features_extracted] = mainFunc_radiomicFeatureExtraction(workingDir, relative_path_to_ROI, outputDir, prefix_of_output_csv, renew_file);
end
