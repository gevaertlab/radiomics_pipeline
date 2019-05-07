function glrlm = computeGrayLevelRunLengthMatrix(I, levels)
nLevels = length(levels);
glrlm = zeros(nLevels, max(size(I)), 13);
I = padIwithNAN(I);
glrlm = compute1D(I, levels, glrlm);
glrlm = compute2D(I, levels, glrlm);
glrlm = compute3D(I, levels, glrlm);

function glrlm = compute3D(I, levels, glrlm)
directions = [1 1 1; -1 1 1; 1 -1 1; -1 -1 1];
directionIndex = 10:13;
nDirections = size(directions, 1);
for i = 1:nDirections
    currentDirection = directions(i,:);
    I = reorientImage(I, currentDirection);
    glrlm(:,:,directionIndex(i)) = computeGLRLM3D(I, levels);
end

function glrlm = computeGLRLM3D(I, levels)
nLevels = length(levels);
nRow = size(I,1);
nCol = size(I,2);
nSli = size(I,3);
glrlm = zeros(nLevels, max(size(I)));
for i = 1:nRow
    for j = 1:nCol
        r = i:(nRow*(i >= j)+(nCol-j+i)*(j > i));
        c = j:((nRow-i+j)*(i >= j)+(nCol*(j > i)));
        s = 1:((nRow-i+1)*(i >= j)+(nCol-j+1)*(j > i));
        voxelIndex = sub2ind(size(I), r, c, s);
        intensityVector = I(voxelIndex);
        intensityVector(isnan(intensityVector)) = [];
        if isempty(intensityVector)
            continue
        else
            index = [find(intensityVector(1:end-1) ~= intensityVector(2:end)), length(intensityVector)];
            len = diff([0 index]);
            val = intensityVector(index);
            if isempty(val)
                continue
            else
                n = size(val);
                idx = zeros(n);
                for k = 1:n(2)
                    idx(k) = find(levels == val(k));
                end
                glrlm = glrlm+accumarray([idx; len]', 1, [nLevels max(size(I))]);
            end
        end
    end
end
for i = 1:nRow
    for j = 2:nSli
        r = i:(nRow*(i >= j)+(nSli-j+i)*(j > i));
        c = 1:((nRow-i+1)*(i >= j)+(nSli-j+1)*(j > i));
        s = j:((nRow-i+j)*(i >= j)+nSli*(j > i));
        voxelIndex = sub2ind(size(I), r, c, s);
        intensityVector = I(voxelIndex);
        intensityVector(isnan(intensityVector)) = [];
        if isempty(intensityVector)
            continue
        else
            index = [find(intensityVector(1:end-1) ~= intensityVector(2:end)), length(intensityVector)];
            len = diff([0 index]);
            val = intensityVector(index);
            if isempty(val)
                continue
            else
                n = size(val);
                idx = zeros(n);
                for k = 1:n(2)
                    idx(k) = find(levels == val(k));
                end
                glrlm = glrlm+accumarray([idx; len]', 1, [nLevels max(size(I))]);
            end
        end
    end
end
for i = 2:nCol
    for j = 2:nSli
        r = 1:((nCol-i+1)*(i >= j)+(nSli-j+1)*(j > i));
        c = i:(nCol*(i >= j)+(nSli-j+i)*(j > i));
        s = j:(nCol-i+j)*(i >= j)+(nSli*(j > i));
        voxelIndex = sub2ind(size(I), r, c, s);
        intensityVector = I(voxelIndex);
        intensityVector(isnan(intensityVector)) = [];
        index = [find(intensityVector(1:end-1) ~= intensityVector(2:end)), length(intensityVector)];
        if isempty(intensityVector)
            continue
        else
            len = diff([0 index]);
            val = intensityVector(index);
            if isempty(val)
                continue
            else
                n = size(val);
                idx = zeros(n);
                for k = 1:n(2)
                    idx(k) = find(levels == val(k));
                end
                glrlm = glrlm+accumarray([idx; len]', 1, [nLevels max(size(I))]);
            end
        end
    end
end

function glrlm = compute2D(I, levels, glrlm)
directions = [1 1 0; -1 1 0; 0 1 1; 0 -1 1; 1 0 1; -1 0 1];
directionIndex = 4:9;
nDirections = size(directions, 1);
for i = 1:nDirections
    currentDirection = directions(i,:);
    I = reorientImage(I, currentDirection);
    glrlm(:,:,directionIndex(i)) = computeGLRLM2D(I, levels);
end

function glrlm = computeGLRLM2D(I, levels)
nRow = size(I,1);
nLevels = length(levels);
glrlm = zeros(nLevels, nRow);
for i = 1:nRow
    II = I(:,:,i);
    for j = (-nRow+1):(nRow-1)
        intensityVector = diag(II, j);
        intensityVector(isnan(intensityVector)) = [];
        if isempty(intensityVector)
            continue
        else
            intensityVector = reshape(intensityVector, [1 length(intensityVector)]);
            index = [find(intensityVector(1:end-1) ~= intensityVector(2:end)), length(intensityVector)];
            len = diff([0 index]);
            val = intensityVector(index);
            if isempty(val)
                continue
            else
                n = size(val);
                idx = zeros(n);
                for k = 1:n(2)
                    idx(k) = find(levels == val(k));
                end
                glrlm = glrlm+accumarray([idx; len]', 1, [nLevels max(size(I))]);
            end
        end
    end
end

function I = reorientImage(I, currentDirection)
if currentDirection(1) == -1
    I = flip(I,1);
end
if currentDirection(2) == -1
    I = flip(I,2);
end
if currentDirection(1) == 0
    I = permute(I, [2 3 1]);
end
if currentDirection(2) == 0
    I = permute(I, [1 3 2]);
end

function glrlm = compute1D(I, levels, glrlm)
% directions = [1 0 0; 0 1 0; 0 0 1];
nRow = size(I,1);
nCol = size(I,2);
nSli = size(I,3);
nLevels = length(levels);

%% 1 0 0
for i = 1:nSli
    for j = 1:nRow
        intensityVector = squeeze(I(j,:,i));
        intensityVector(isnan(intensityVector)) = [];
        if isempty(intensityVector)
            continue
        else
            intensityVector = reshape(intensityVector, [1 length(intensityVector)]);
            index = [find(intensityVector(1:end-1) ~= intensityVector(2:end)), length(intensityVector)];
            len = diff([0 index]);
            val = intensityVector(index);
            if isempty(val)
                continue
            else
                n = size(val);
                idx = zeros(n);
                for k = 1:n(2)
                    idx(k) = find(levels == val(k));
                end
                glrlm(:,:,1) = glrlm(:,:,1)+accumarray([idx; len]', 1, [nLevels max(size(I))]);
            end
        end
    end
end

%% 0 1 0
for i = 1:nSli
    for j = 1:nCol
        intensityVector = squeeze(I(:,j,i));
        intensityVector(isnan(intensityVector)) = [];
        if isempty(intensityVector)
            continue
        else
            intensityVector = reshape(intensityVector, [1 length(intensityVector)]);
            index = [find(intensityVector(1:end-1) ~= intensityVector(2:end)), length(intensityVector)];
            len = diff([0 index]);
            val = intensityVector(index);
            if isempty(val)
                continue
            else
                n = size(val);
                idx = zeros(n);
                for k = 1:n(2)
                    idx(k) = find(levels == val(k));
                end
                glrlm(:,:,2) = glrlm(:,:,2)+accumarray([idx; len]', 1, [nLevels max(size(I))]);
            end
        end
    end
end

%% 0 0 1
for i = 1:nRow
    for j = 1:nCol
        intensityVector = squeeze(I(i,j,:));
        intensityVector(isnan(intensityVector)) = [];
        if isempty(intensityVector)
            continue
        else
            intensityVector = reshape(intensityVector, [1 length(intensityVector)]);
            index = [find(intensityVector(1:end-1) ~= intensityVector(2:end)), length(intensityVector)];
            len = diff([0 index]);
            val = intensityVector(index);
            if isempty(val)
                continue
            else
                n = size(val);
                idx = zeros(n);
                for k = 1:n(2)
                    idx(k) = find(levels == val(k));
                end
                glrlm(:,:,3) = glrlm(:,:,3)+accumarray([idx; len]', 1, [nLevels max(size(I))]);
            end
        end
    end
end

function II = padIwithNAN(I)
II = nan(max(size(I)),max(size(I)),max(size(I)));
II(1:size(I,1),1:size(I,2),1:size(I,3)) = I;