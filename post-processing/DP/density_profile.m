%% Water Density Profile (Above/Below Graphene)
%
% F. Tarulli – Politecnico di Torino, Mar 20, 2025
%
% Computes the z-resolved water density (g/cm³) from atom counts in
% `counters.mat`.

%%
clc;
clear vars
tic;

%% 
%Layers
maxZheight = 15;                    % MUST be consistent with the same variable in 'data_parser.m' script
thickness = 0.1;                    % MUST be consistent with the same variable in 'data_parser.m' script
maxLayers = maxZheight/thickness; 


    filename = '../../lammps/DP/water-graph_density.dump';
    fileID = fopen(filename, 'r');
    try
        for i=1:5
            fgets(fileID);
        end
        line = textscan(fileID, '%f %f' ,1);
        line = cell2mat(line);
        widthBOX = line(end) - line(1);
    
        line = textscan(fileID, '%f %f' ,1);
        line = cell2mat(line);
        lengthBOX = line(end) - line(1);
        fclose(fileID);
    catch
        ME = MException('Warn:AreaNotFound','Cross sectional area not found. File .dump missing: import "widthBOX" and "lengthBOX" manually\n');
        throw(ME)
    end


N_A = 6.022140857e23;
MM_O = 15.999;
MM_H = 1.0079;

V = widthBOX * lengthBOX * thickness * 1e-24;        % Approximatated bins volume

atom_types = struct('name', {'O', 'H',}, 'index', {1, 2});
molar_masses = struct('O', MM_O, 'H', MM_H);


counter = load('Counters/counters.mat').counterDOWN;

mean_counters = struct();
masses = struct();
sumasses=zeros(maxLayers,1);

for i = 1:length(atom_types)
    atom_name = atom_types(i).name;
    counter_name = [atom_name '_count'];
    
    mean_counters.(atom_name) = mean(counter.(counter_name), 2);
    
    masses.(atom_name) = mean_counters.(atom_name) / N_A * molar_masses.(atom_name);
    sumasses = sumasses + masses.(atom_name);
end

densityDOWN=sumasses/V;


counter = load('Counters/counters.mat').counterUP;

mean_counters = struct();
masses = struct();
sumasses=zeros(maxLayers,1);

for i = 1:length(atom_types)
    atom_name = atom_types(i).name;
    counter_name = [atom_name '_count'];
    
    mean_counters.(atom_name) = mean(counter.(counter_name), 2);

    masses.(atom_name) = mean_counters.(atom_name) / N_A * molar_masses.(atom_name);
    sumasses = sumasses + masses.(atom_name);
end

densityUP=sumasses/V;

density = [flip(densityDOWN);densityUP];

%%
fprintf('Elapsed time: %.1f s\n', toc);

%% Plot

position = (thickness/2:thickness:maxZheight)';
positions = [-flip(position);position];

x_fine = linspace(positions(1), positions(end), 500); 
y_fine = interp1(positions, density, x_fine, 'pchip'); 

figure;
plot(x_fine, y_fine, 'r', 'LineWidth', 2);
hold on
plot(linspace(-maxZheight,maxZheight, length(density)), density, 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
title(sprintf('Water Density profile'))
xlabel('Z (Å)')
ylabel('Density (g/cm^3)')


%% (Optionally) save to plot in python
%save([path, '/', num2str(oxid) '%.mat'],  'x_fine', 'positions', 'y_fine', 'density');
