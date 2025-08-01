% water density with Molecular Surface generated with EDTSurf 
% Last code (8/7/25) for Linux

%%
clear vars; close all; clc
tic

%% Settings
oxid = ;                                      % Oxidation degree

path_dump = fullfile('../../lammps/DP');
path_ply = fullfile('MS','ply_ms');

atom_types = struct('name', {'O','H'}, 'index', {1,2});

N = 78148 + 2 * round(oxid/100 * 1972);       % 78148 water particles | 1972 graphene carbon atoms

warning('off') % ('Warning: Duplicate data points have been detected and removed - corresponding values have been averaged.' during interpolation)

% Simulation details
acquisition_steps  = 1000;                    % Trajectory (dump) file is saved every this many simulation steps
initial_timestep   = 1;
final_timestep     = 1501;                    % Final frame index; density profile will be averaged from initial_timestep to this frame
thickness          = 0.1;                     % Thickness of layers
maxZheight         = 15;                      % Maximum distance (in Angstrom) above and below the graphene for evaluating the density profile
maxLayers          = maxZheight/thickness;


%% Parsing
dump_file = fullfile(path_dump, 'water-graph_density.dump');
fileID    = fopen(dump_file, 'r');

s = 1; 
for i = 1:length(atom_types)
    atom_name = atom_types(i).name;
    counterUP.([atom_name '_count'])   = zeros(maxLayers, final_timestep - initial_timestep + 1);
    counterDOWN.([atom_name '_count']) = zeros(maxLayers, final_timestep - initial_timestep + 1);
end
V = zeros(final_timestep - initial_timestep + 1, 1);
fprintf('Timestep:\n')

for k = initial_timestep:final_timestep
    timestep = k - 1;
    fprintf('%d\n',timestep)        

    if k==initial_timestep
        for i=1:9*k+timestep*N 
            fgets(fileID);
        end
    end
        dati = textscan(fileID, '%f %f %f %f %f', N);
    for i=1:10
        fgets(fileID);
    end

    dati = cell2mat(dati);
    if size(dati,1)~=N
        fprintf('ERROR! timestep=%d\n',timestep)
    end

%% Processing

    C_data = dati(dati(:, 2) == 3 | dati(:, 2) == 4, 3:5);   % (type=3,4)
    H_data = dati(dati(:, 2) == 6,             3:5);         % (type=6)

    F_C = scatteredInterpolant(C_data(:,1), C_data(:,2), C_data(:,3), 'nearest');
    F_H = scatteredInterpolant(H_data(:,1), H_data(:,2), H_data(:,3), 'nearest');

    % Reading mesh PLY 
    file_path_full = fullfile(path_ply,sprintf('ms_%d.ply', timestep));
    if ~exist(file_path_full,'file')
        fprintf('WARNING: %s do not exist!\n', file_path_full );
        continue
    end
    [faces, verts] = fastPlyRead(file_path_full);
    zC_vert = F_C(verts(:,1), verts(:,2)); 
    dist    = verts(:,3) - zC_vert;

    %-----------------------------
    %   Above graphene
    aboveIdx    = (dist > 0);                   
    newVertsIdx = find(aboveIdx);               
    mapOldToNew = zeros(size(aboveIdx));
    mapOldToNew(newVertsIdx) = 1:numel(newVertsIdx);

    faceMask       = aboveIdx(faces(:,1)) & aboveIdx(faces(:,2)) & aboveIdx(faces(:,3));
    facesFiltered  = faces(faceMask, :);
    facesFiltered  = [mapOldToNew(facesFiltered(:,1)), ...
                      mapOldToNew(facesFiltered(:,2)), ...
                      mapOldToNew(facesFiltered(:,3))];
    vertsFiltered  = verts(newVertsIdx, :);
    F_ply_up       = scatteredInterpolant(vertsFiltered(:,1), ...
                                          vertsFiltered(:,2), ...
                                          vertsFiltered(:,3), 'nearest','none');
    %-----------------------------

    %-----------------------------
    %   below graphene
    belowIdx           = (dist < 0);
    newVertsIdx_down   = find(belowIdx);
    mapOldToNew_down   = zeros(size(belowIdx));
    mapOldToNew_down(newVertsIdx_down) = 1:numel(newVertsIdx_down);

    faceMask_down      = belowIdx(faces(:,1)) & belowIdx(faces(:,2)) & belowIdx(faces(:,3));
    facesFiltered_down = faces(faceMask_down, :);
    facesFiltered_down = [mapOldToNew_down(facesFiltered_down(:,1)), ...
                          mapOldToNew_down(facesFiltered_down(:,2)), ...
                          mapOldToNew_down(facesFiltered_down(:,3))];
    vertsFiltered_down = verts(newVertsIdx_down, :);
    F_ply_down         = scatteredInterpolant(vertsFiltered_down(:,1), ...
                                              vertsFiltered_down(:,2), ...
                                              vertsFiltered_down(:,3), 'nearest','none');
    %-----------------------------

    z_vals_UP   = F_ply_up(dati(:,3), dati(:,4));
    z_vals_DOWN = F_ply_down(dati(:,3), dati(:,4));

    typeAtom = dati(:,2);
    zCoord   = dati(:,5);

    %-----------------------------
    %   Above graphene
    zUP           = zCoord - z_vals_UP;
    layerIndexUp  = floor(zUP / thickness) + 1;  
    validUp       = layerIndexUp >= 1 & layerIndexUp <= maxLayers; 
    %-----------------------------

    %-----------------------------
    %   below graphene
    zDOWN            = z_vals_DOWN - zCoord;
    layerIndexDown   = floor(zDOWN / thickness) + 1;
    validDown        = layerIndexDown >= 1 & layerIndexDown <= maxLayers;
    %-----------------------------


    for iType = 1:length(atom_types)
        atom_name   = atom_types(iType).name;
        atom_number = atom_types(iType).index;

        maskAtom = (typeAtom == atom_number);

        maskUp              = validUp & maskAtom;
        idxLayersUp         = layerIndexUp(maskUp); 
        countUP             = accumarray(idxLayersUp, 1, [maxLayers, 1]);
        counterUP.([atom_name '_count'])(:, s) = countUP;

        maskDown            = validDown & maskAtom;
        idxLayersDown       = layerIndexDown(maskDown);
        countDOWN           = accumarray(idxLayersDown, 1, [maxLayers, 1]);
        counterDOWN.([atom_name '_count'])(:, s) = countDOWN;
    end


    s = s + 1;
end
fclose(fileID);

%% Saving 
% number of water molecules in each layer for each timestep
mkdir('counters')
save(fullfile('counters', 'counters.mat'), 'counterUP', 'counterDOWN');

fprintf('Elapsed time= %f',toc)
