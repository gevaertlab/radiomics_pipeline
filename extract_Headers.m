function [features_extracted] = extract_Headers(workingDir, relative_path_to_ROI, outputDir, prefix_of_output_csv, renew_file)
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

%% collect dicom headers for QC, mask and volume for later use. 
% Slice spacing is computed according to: <https://github.com/mvallieres/radiomics/>
dicomInfo = struct;
% column 1, all slices' headers; column 2, computed sliceS(slice spacing).
% column 3, slice thickness.
for j =1:length(All_Data_Path) %
    fprintf('Reading headers, mask and volume from %d th case \n',j);
    %load data
    Directory = fullfile(filepath, All_Data_Path(j).name);
    elements=subdir(Directory); % list folders and subfolders. subdir()
    % Copyright (c) 2015 Kelly Kearney.
    nElements = length(elements); % includes all elements:e.g. .DS_Store files, mask file.
    volume = cell(1,1,nElements); % initialize the volume depth the same as the number of all elements. will trim later.
    dicomHeaders = []; % initialization, dicomHeaders refer to headers for
    % all DICOM files from a series
    SEG = []; % will store mask array from DSO file. The modality is "SEG".
    
    dicomInfo(j).patID = All_Data_Path(j).name; % patient ID
    
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
%                     fprintf('\nID:%s, sliceNumber:%s shape:%s \n',dicomInfo(j).patID, num2str(sliceNumber), num2str(shape(1)));
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
        
        dicomInfo(j).SEGvol = SEGvol;
        dicomInfo(j).SEGheader = SEGheader;
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
    dicomInfo(j).volume = volume;
    dicomInfo(j).SEGvol = SEGvol;
    
    dicomInfo(j).headers = dicomHeaders;
    dicomInfo(j).dicomInfo1slice = dicomHeaders(1);
    dicomInfo(j).SliceSpacing = SliceSpacing;
    dicomInfo(j).SliceThickness = dicomHeaders(1).SliceThickness;
    if isfield(dicomHeaders,'RescaleType')
        dicomInfo(j).RescaleType = dicomHeaders(1).RescaleType; % not all DICOM
    % series has this fieldname. TCGA collection "TCGA-HNSC" does, while
    % "NSCLC-Radiogenomics" does not.
    end
    dicomInfo(j).RescaleSlope = dicomHeaders(1).RescaleSlope;
    dicomInfo(j).RescaleIntercept = dicomHeaders(1).RescaleIntercept;
end % end reading all patients
% dicomData will contain data from 2 patients.

%% check if the patIDs are unique.
patIDs = {dicomInfo.patID};
if length(unique(patIDs)) == size(dicomInfo,2)
    fprintf('All patIDs are unique, ready to move on\n');
else
    fprintf('One or more patIDs are duplicated, PLEASE correct the ID\n');
end
% end
end