function [structureArray] = appendStruct(structureArray,newStructure)

if isempty(structureArray)
    structureArray = newStructure;
    return
end

structLength = length(structureArray);
fields = fieldnames(structureArray(1));
nFields = length(fields);

for i = 1:nFields
    try
        structureArray(structLength + 1).(fields{i}) = newStructure.(fields{i});
    catch
        structureArray(structLength + 1).(fields{i}) = 'FIELD NOT PRESENT';
    end
end

end