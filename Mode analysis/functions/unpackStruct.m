function unpack(inputStruct)
    % Check if the input is a structure
    if ~isstruct(inputStruct)
        error('Input must be a structure.');
    end

    % Get the field names of the structure
    fields = fieldnames(inputStruct);

    % Assign each field to a variable in the caller workspace
    for i = 1:length(fields)
        fieldName = fields{i};
        fieldValue = inputStruct.(fieldName);
        assignin('caller', fieldName, fieldValue);
    end
end