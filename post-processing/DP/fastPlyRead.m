function [faces, verts] = fastPlyRead(filename)
% fastPlyRead - Fast reader for ASCII PLY files with vertex colors.
%
%   [faces, verts] = fastPlyRead(filename) reads the PLY file specified by
%   filename and returns:
%       verts - an n-by-3 matrix of vertex coordinates (x, y, z)
%       faces - an m-by-3 matrix with indices of vertices forming each face.
%
% The function expects the PLY file to have:
%   - A header that specifies an element "vertex" with 6 properties:
%         x, y, z, red, green, blue.
%   - An element "face" where each line begins with the number of vertices 
%     (assumed to be 3), followed by three vertex indices, and then three
%     color values (red, green, blue) which are ignored.
%
% Face indices in the file are assumed to be 0-indexed and are converted to
% MATLAB’s 1-indexed format.

fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open the file %s', filename);
end

numVerts = 0;
numFaces = 0;
vertexPropCount = 0;
inVertexSection = false;

% --- Read header ---
while true
    line = fgetl(fid);
    if line == -1
        error('Unexpected end of file while reading header.');
    end
    
    if contains(line, 'element vertex')
        numVerts = sscanf(line, 'element vertex %d');
        inVertexSection = true;
    elseif contains(line, 'element face')
        numFaces = sscanf(line, 'element face %d');
        inVertexSection = false;  % finished vertex section
    elseif startsWith(strtrim(line), 'property') && inVertexSection
        vertexPropCount = vertexPropCount + 1;
    elseif strcmp(strtrim(line), 'end_header')
        break;
    end
end

% --- Read vertex data ---
% Build a format string based on the number of vertex properties.
% For your file, vertexPropCount should be 6.
vertexFormat = repmat('%f ', 1, vertexPropCount);
% Read all vertex data in one call.
vertsData = textscan(fid, vertexFormat, numVerts);
vertsData = cell2mat(vertsData);
% Extract the first three columns (x, y, z).
verts = vertsData(:, 1:3);

% --- Read face data ---
% Each face line is structured as:
%    vertex_count v1 v2 v3 red green blue
% We ignore vertex_count and the face color by reading:
faceFormat = '%*d %d %d %d %*d %*d %*d';
facesData = textscan(fid, faceFormat, numFaces);
faces = cell2mat(facesData);

% Convert face indices from 0-indexed to 1-indexed.
if ~isempty(faces) && min(faces(:)) == 0
    faces = faces + 1;
end

fclose(fid);
end
