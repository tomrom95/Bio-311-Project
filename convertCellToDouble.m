function [out] = convertCellToDouble( cellTable )
    for i=1:length(cellTable)
        for j=1:length(cellTable(1,:))
            if iscellstr(cellTable(i,j))
                cellTable{i,j} = str2double(cellTable(i,j));
            end
        end
    end
    out = cellTable;
end

