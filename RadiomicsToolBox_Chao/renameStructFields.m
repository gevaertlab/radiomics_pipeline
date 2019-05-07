function structName = renameStructFields(structName, type, name)
% type - like 'HHH', 'HHL' of wavelet decompositions.
% name - like 'wavelet'.
oldNames = fieldnames(structName);
nFields = length(oldNames);
for i = 1 :nFields
    newName = [name type oldNames{i}];
    structName.(newName) = structName.(oldNames{i}); % if "Scalar structure required for this assignment", switch to [structName.(newName)] = structName.(oldNames{i})
    structName = rmfield(structName, oldNames{i});
end