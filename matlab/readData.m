function [data, geneList, cell_types, cell_IDs, cell_IDs_unique] = readData(fileDir)

str = fileread(fileDir);   %read entire file into string
parts = strtrim(regexp( str, '(\r|\n)+', 'split'));  %split by each line
cell_IDs = strtrim( regexp(parts{1}, '\s+', 'split'));  %columns
% cell_IDs(1) = [];
ncol = length(cell_IDs);  %number of columns
parts(1)= [];  %remove column headers
nrows = length(parts);  %number of rows
geneList = cell(nrows-1, 1);  %pre-allocate empty cell array for data
cellCluster = cell(nrows-1, 1);  %pre-allocate empty cell array for data
data = zeros(nrows-1, ncol);

for k=1:nrows-1
    buffer = strtrim(regexp( parts{k}, '\s+', 'split'));   %split by spaces
    buffer = strrep(buffer, '"', '');
    geneList{k} = buffer{1};
    C = buffer(2:end);
    S = sprintf('%s ', C{:});
    data(k,:) = sscanf(S, '%f');
end
geneList = strrep(geneList, '"', '');
cell_IDs_unique = unique(cell_IDs);
cell_types = ones(size(cell_IDs));
for i = 1:length(cell_IDs_unique)
    idx = find(~cellfun(@isempty,strfind(cell_IDs,cell_IDs_unique{i})));
    cell_types(idx) = i;
end
