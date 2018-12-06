NCols = 9;
NRows = 35299;
file = 'w:\Bladed\2Benergy_MC141\powprod_DOFdeactivated_binary.$04';
fid = fopen(file, 'r');
data = fread(fid,NRows*NCols,'float32');
data = (reshape(data, NCols, NRows))';

