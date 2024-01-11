% Reads in a list of points in a text file with the format used
% in Camino.
% author: Daniel C. Alexander (d.alexander@ucl.ac.uk)
% $Id $


function pts = ReadCaminoElecPS(pointsFile)

fid = fopen(pointsFile, 'r', 'b');
A = fscanf(fid, '%f');
fclose(fid);

% The first number is the number of points in the pointset.
noPoints = A(1);
A = A(2:end);

% Get the points in a sensible data structure.
pts = reshape(A, 3, noPoints)';


