% extract features


%% TCGA-TCIA data
workingDir = '/Users/messi/Documents/private_research/Stanford_BMI/HeadNeck_Cancer/data/original/hnscc_CT/';
cd (workingDir);

relative_path_to_ROI = 'hnscc_113pats_ROI3_20171213_byMurilo/Patients_DICOM_ROI_2/';
prefix_of_output_csv = 'TCGA_TCIA_radiomicFeatures_'; % date will auto added as suffix of the file name.
outputDir = workingDir;
renew_file = '/Users/messi/Documents/private_research/Stanford_BMI/HeadNeck_Cancer/data/original/hnscc_CT/TCGA_TCIA_radiomicFeatures_20190504.mat'; % file to be appended with feature rows for new patients
[features_extracted] = mainFunc_radiomicFeatureExtraction(workingDir, relative_path_to_ROI, outputDir, prefix_of_output_csv, renew_file);

%% Stanford data
workingDir = '/Users/messi/Documents/private_research/Stanford_BMI/HeadNeck_Cancer/data/original/hnscc_CT/';
cd (workingDir);

relative_path_to_ROI = 'Stanford_Cohort_HNSC_from_Murilo_20180328/ROI_DI/';
prefix_of_output_csv = 'Stanford_radiomicFeatures_';
outputDir = workingDir;
renew_file = '';
[features_extracted] = mainFunc_radiomicFeatureExtraction(workingDir, relative_path_to_ROI, outputDir, prefix_of_output_csv, renew_file);
