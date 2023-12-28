function Matrix = ReadBinMat(Filename)
% This function is to read double type matrix stored in binary file
fid = fopen(Filename, 'rb');

% read rows and columns
rows = fread(fid, 1, 'int');
cols = fread(fid, 1, 'int');

% read matrix data
Matrix = fread(fid, [rows, cols], 'double');

fclose(fid);

end

