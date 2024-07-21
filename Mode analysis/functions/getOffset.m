function offset = getOffset(sublattice_size)
    size_of_rows = [0:sublattice_size-2];
    for i = 1:sublattice_size
        size_of_rows = [size_of_rows, sublattice_size-1];
    end
    offset = size_of_rows;
end
