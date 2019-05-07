function [features_extracted] = mainFunc_radiomicFeatureExtraction(workingDir, relative_path_to_ROI, outputDir, prefix_of_output_csv, renew_file)
% Author: Chao Huang (huangchao312@gmail.com, huangchao09@zju.edu.cn)
% updated: April 4th, 2019.
% Sourcefile: mainFile_feature_extraction_from_DICOM_n_DSO.m

addpath(genpath('/Users/sarahmattonen/Documents/Projects/QIFP/ChaosPipeline/radiomicFeatureExtractionPipeline/'))
% Inputs:
%   - workingDir:  quoted in ''.the working directory, where you store the images, DSO files and
%   want to store the results.
%   - relative_path_to_ROI:  quoted in ''.the path to ROI folder relative to the
%   workingDir. This folder only have all the images for the studied
%   patients, with images of each patient stored in a separate subfolder.
%   - outputDir:  quoted in ''.
%   - prefix_of_output_csv: quoted in ''.
%   - renew_file:  quoted in ''. Please pass '' to do from scratch. Otherwise, if the execution breaks at some point. you can
%   renew the process by passing the newest '.mat' file path.

% workingDir =
% '/Users/messi/Documents/private_research/Stanford_BMI/TCIA_TCGA/';
cd (workingDir);

% relative_path_to_ROI = 'NSCLC Radiogenomics/';
% prefix_of_output_csv = 'NSCLC_radiogenomics_radiomicFeatures_';


%% initializers 

filepath = fullfile(workingDir, relative_path_to_ROI);
All_Data_Path = dir(filepath);
% delete '.', '..', '.DS_Store', and '_DS_Store'
toRemoveIdx = [];
for i = 1:length(All_Data_Path)
    if startsWith(All_Data_Path(i).name,'.') || startsWith(All_Data_Path(i).name,'..') || strcmp(All_Data_Path(i).name,'.DS_Store') || strcmp(All_Data_Path(i).name,'_DS_Store') || strcmp(All_Data_Path(i).name,'__.DS_Store') || strcmp(All_Data_Path(i).name,'._.DS_Store') || contains(All_Data_Path(i).name,'.csv') || contains(All_Data_Path(i).name,'.xml') || contains(All_Data_Path(i).name,'.mat') % '._.DS_Store' required for Linux.
        toRemoveIdx = [toRemoveIdx,i];
    end;
end;
All_Data_Path(toRemoveIdx,:) = [];

formatOut = 'yyyymmdd'; % format date time to be like '20171225'

% end

%% used for debugging
% All_Data_Path = All_Data_Path([1],:);

%% start to extract radiomic features
if isempty(renew_file)
    radiomicFeatures=[];
    starter = 1;
else
    load(renew_file, 'radiomicFeatures');
    starter = length(radiomicFeatures)+1;
end
% load('TCGA_TCIA_radiomicFeatures_20181106.mat', 'radiomicFeatures');

for j = starter:length(All_Data_Path) %
    % collect dicom headers for this case
    dicomInfo = struct;
        %load data
    Directory = fullfile(filepath, All_Data_Path(j).name);
    elements=subdir(Directory); % list folders and subfolders. subdir()
    % Copyright (c) 2015 Kelly Kearney.
    nElements = length(elements); % includes all elements:e.g. .DS_Store files, mask file.
    volume = cell(1,1,nElements); % initialize the volume depth the same as the number of all elements. will trim later.
    dicomHeaders = []; % initialization, dicomHeaders refer to headers for
    % all DICOM files from a series
    SEG = []; % will store mask array from DSO file. The modality is "SEG".
    
    dicomInfo.patID = All_Data_Path(j).name; % patient ID
    
    % starts to read
    sliceNumber = 0;
    for elementNumber = 1:nElements
        elementName = elements(elementNumber).name;
        [p,name,ext] = fileparts(elementName);
        baseName=strcat(name,ext);
        if strcmp(ext, '.dcm') && ~startsWith(baseName,'.') && ~startsWith(baseName,'..') && ~contains(baseName,'.DS_Store') && ~contains(baseName,'_DS_Store') && ~contains(baseName,'__.DS_Store') && ~contains(baseName,'.csv') && ~contains(baseName,'.xml') && ~contains(baseName,'.mat')
            % Good enough for Linux, add conditions for MAC and Windows.
            % strcmp() compares indentical strings. contains() means matching a
            % pattern.
            elementFullFile = elementName;
            if isdicom(elementFullFile)
                tmp = dicominfo(elementFullFile);
                if strcmp(tmp.Modality,'RTSTRUCT')
                    RTstruct = tmp;
                elseif strcmp(tmp.Modality,'REG')
                    REG = tmp;
                elseif strcmp(tmp.Modality,'SEG')
                    SEGheader = tmp;
                    SEG = squeeze(dicomread(elementFullFile));
                elseif strcmp(tmp.Modality,'MR') || strcmp(tmp.Modality,'PT') || strcmp(tmp.Modality,'CT')
                    sliceNumber = sliceNumber + 1;
                    volume{sliceNumber} = int16(dicomread(elementFullFile)); %double(dicomread(elementFullFile))
%                     shape = size(dicomread(elementFullFile));
%                     fprintf('\nID:%s, sliceNumber:%s shape:%s \n',dicomInfo.patID, num2str(sliceNumber), num2str(shape(1)));
                    dicomHeaders = appendStruct(dicomHeaders,tmp); % appendStruct()
                    % refer <https://github.com/mvallieres/radiomics/>.
                end
            end
        end
    end % end reading for one patient
    nSlices = sliceNumber; % Total number of slices
    volume = volume(1:nSlices); % Suppress empty cells in images
    
    % DETERMINE THE SCAN ORIENTATION
    dist = [abs(dicomHeaders(2).ImagePositionPatient(1) - dicomHeaders(1).ImagePositionPatient(1)), ...
        abs(dicomHeaders(2).ImagePositionPatient(2) - dicomHeaders(1).ImagePositionPatient(2)), ...
        abs(dicomHeaders(2).ImagePositionPatient(3) - dicomHeaders(1).ImagePositionPatient(3))];
    [~,index] = max(dist);
    if index == 1
        orientation = 'Sagittal';
    elseif index == 2
        orientation = 'Coronal';
    else
        orientation = 'Axial';
    end
    
    % For original DICOM image series, SORT THE IMAGES AND DICOM HEADERS after
    % extracting the right order indices
    slicePositions = zeros(1,nSlices);
    for sliceNumber = 1:nSlices
        slicePositions(sliceNumber) = dicomHeaders(sliceNumber).ImagePositionPatient(index);
    end
    [~,indices_img] = sort(slicePositions);
    
    % transform values stored as uint16 to int16, and to be HU(Housefield Units)
    volume = int16(cell2mat(volume(indices_img)));
    type = dicomHeaders(1).Modality;
    if strcmp(type,'PT') || strcmp(type,'CT')
        if strcmp(type,'PT')
            type = 'PET';
        end
        for i=1:size(volume,3)
            volume(:,:,i)=volume(:,:,i)*dicomHeaders(i).RescaleSlope + dicomHeaders(i).RescaleIntercept;
        end
    end

    dicomHeaders = dicomHeaders(indices_img);
    
    if exist('SEGheader', 'var') && isstruct(SEGheader)
        slicePositions_DSO = zeros(1,nSlices);
        for nSDSO = 1:nSlices
            tmpSDSO = SEGheader.ReferencedSeriesSequence. ...
                Item_1.ReferencedInstanceSequence.(['Item_' num2str(nSDSO)]). ...
                ReferencedSOPInstanceUID;
            
            slice_index = structfind(dicomHeaders,'SOPInstanceUID',tmpSDSO); % find the slice index in dicomHeaders to identify the img slice corresponding to this DSO slice
            slicePositions_DSO(nSDSO) = dicomHeaders(slice_index).ImagePositionPatient(index);
        end
        [~,indices_DSO] = sort(slicePositions_DSO);
        
        % NOTE: slicePositions and slicePositions_DSO should have the same values
        % set, but probably have different order of the values. That's why we need
        % to sort them separately to get the dicom images volume and the mask
        % volume with the same slice orders.
        
        SEGvol = SEG(:,:,indices_DSO);
        
        dicomInfo.SEGvol = SEGvol;
        dicomInfo.SEGheader = SEGheader;
    else
        % sometimes only ROI are given where pixels outside the ROI are set
        % to 0 while pixels inside the ROI are the original intensities of
        % the images.
        % 
        SEGvol = volume;
        SEGvol(SEGvol~=0)=1;
    end
    

    % calculate slice spacing
    s1 = round(0.5*nSlices); s2 = round(0.5*nSlices) + 1; % Slices selected to
    % calculate slice spacing
    SliceSpacing =  sqrt(sum((dicomHeaders(s1).ImagePositionPatient - dicomHeaders(s2).ImagePositionPatient).^2));
    dicomInfo.volume = volume;
    dicomInfo.SEGvol = SEGvol;
    
    dicomInfo.headers = dicomHeaders;
    dicomInfo.dicomInfo1slice = dicomHeaders(1);
    dicomInfo.SliceSpacing = SliceSpacing;
    dicomInfo.SliceThickness = dicomHeaders(1).SliceThickness;
    if isfield(dicomHeaders,'RescaleType')
        dicomInfo.RescaleType = dicomHeaders(1).RescaleType; % not all DICOM
    % series has this fieldname. TCGA collection "TCGA-HNSC" does, while
    % "NSCLC-Radiogenomics" does not.
    end
    dicomInfo.RescaleSlope = dicomHeaders(1).RescaleSlope;
    dicomInfo.RescaleIntercept = dicomHeaders(1).RescaleIntercept;
    
    ID = dicomInfo.headers.PatientID;
    fprintf('\n--------------Extracting features for %d th case(%s)...\n', j, ID);
    DataVolume = dicomInfo.volume; 
    MaskVolume = dicomInfo.SEGvol;
  
    %feature extraction
    %--------------------------------------------------------------------------
    % These parameters are mostly for prepareVolume(volume,mask,scanType,
    % pixelW,sliceS,R,scale,textType,quantAlgo,Ng) from 
    % <https://github.com/mvallieres/radiomics/>
    % The values assigments were modified by Chao according to the original
    % DICOM info embedded in the DICOM files.
    Para.type = 'CT'; % scanType: equals to prepareVolume()'s 'scanType': String 
    % specifying the type of scan analyzed. Either 'PETscan', 'MRscan' or 
    % 'Other'.
    Para.pixelW = dicomInfo.headers(1).PixelSpacing(1);% pixelW: pixel width/spacing, Numerical value specifying the 
    % in-plane resolution (mm) of 'volume'.
    Para.sliceS = dicomInfo.SliceSpacing; % sliceS: Numerical value specifying the slice spacing (mm) of 'volume'.
    % Slice Thickness differs from Slice Spacing: more see https://stackoverflow.com/questions/14930222/how-to-calculate-space-between-dicom-slices-for-mpr
    Para.R = 1;% R: Numerical value specifying the ratio of weight to 
    % band-pass coefficients over the weight of the rest of coefficients 
    %(HHH and LLL). Provide R=1 to not perform wavelet band-pass filtering. 
    Para.scale = 1;% scale: Numerical value specifying the scale at which 
    % 'volume' is isotropically resampled (mm). If a string 'pixelW' is 
    % entered as input, the volume will be isotropically resampled at the 
    % initial in-plane resolution of 'volume' specified by 'pixelW'. If 
    % scale_cell = 1, is isotropically resampled to a voxel size of 1x1x1
    % mm3.
    Para.textType = ''; % textType: String specifying for which type of textures 
    % 'volume' is being prepared. Either 'Global' or 'Matrix'. If 'Global',
    % the volume will be prepared for Global texture features computation. 
    % If 'Matrix',the volume will be prepared for matrix-based texture 
    % features computation (i.e. GLCM, GLRLM, GLSZM, NGTDM).
    % NOTE: textType will not be used later. textType will specifically assigned
    % with a string as 'Global' or 'Matrix' to prepareVolume() if
    % necessary.
    Para.quantAlgo = 'Uniform';% quantAlgo: String specifying the
    % quantization algorithm to use on 'volume'. Either 'Equal' for 
    % equal-probability quantization, 'Lloyd' for Lloyd-Max quantization, 
    % or 'Uniform' for uniform quantization. Use only if textType is set to
    % 'Matrix'.
    Para.Ng = 32; % Ng: Integer specifying the number of gray levels in 
    % the quantization process. Use only if textType is set to 'Matrix'.
    
    radiomicFeatures_tmp = radiomicFeatureExtraction_Pipeline(DataVolume,MaskVolume,Para);
    radiomicFeatures_tmp.patID = dicomInfo.headers.PatientID;
    
    radiomicFeatures = [radiomicFeatures,radiomicFeatures_tmp];
    
    SaveFilePath = fullfile(outputDir, [prefix_of_output_csv datestr(now,formatOut) '.mat']);
    save(SaveFilePath,'radiomicFeatures');
    % transform struct to table
    radiomicFeaturesTable = struct2table(radiomicFeatures);

    % optional: move the last column (patID) to the first as conventional practice
    radiomicFeaturesTable = [radiomicFeaturesTable(:,size(radiomicFeaturesTable,2)) radiomicFeaturesTable(:, 1:(size(radiomicFeaturesTable,2)-1))];
    
    SaveTablePath = fullfile(outputDir, [prefix_of_output_csv datestr(now,formatOut) '.csv']);
    writetable(radiomicFeaturesTable, SaveTablePath);

    features_extracted = radiomicFeaturesTable;
    
    fprintf('Done the %d th case, ID:%s \n\n',j,radiomicFeatures_tmp.patID);   
end
% end
end