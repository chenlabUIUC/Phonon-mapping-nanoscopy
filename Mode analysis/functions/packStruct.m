function outputStruct = packStruct(varNames)
    % Check if the input variable names are provided in a cell array
    if ~iscell(varNames)
        error('Variable names must be provided in a cell array.');
    end
    
    % Initialize the output structure
    outputStruct = struct();
    
    % Iterate over the list of variable names
    for i = 1:length(varNames)
        varName = varNames{i};
        
        % Check if the variable exists in the caller workspace
        if evalin('caller', ['exist(''', varName, ''', ''var'')'])
            % Get the value of the variable from the caller workspace
            varValue = evalin('caller', varName);
            % Assign the value to the corresponding field in the structure
            outputStruct.(varName) = varValue;
        else
            warning('Variable "%s" does not exist in the caller workspace.', varName);
        end
    end
end
