%% Graphene OH‐Functionalization
%
% Created by F. Tarulli, Politecnico di Torino, Italy - October 16, 2023
%
% This script reads an existing graphene file (“graphene.lt”), 
% identifies a user‐specified fraction of carbon atoms to functionalize with 
% hydroxyl (–OH) groups, and writes out a new “GrapheneFunct.lt” file 
% incorporating the added O and H atoms, bonds, and angles.
% 
% Workflow:
%    1. Set processing parameters:
%         • dir_path: folder containing input/output files 
%         • oxid: percentage of C atoms to functionalize (0–100)
%         • qCC, qO2, qH2: partial charges for C–OH, O, and H atoms
%         • N: total number of C atoms in the graphene sheet
%    2. Copy the original file and locate the “Data Atoms” section.
%    3. Parse the atom coordinates, types, and charges.
%    4. Compute neighbor distances (k, l) and generate the graphene bond list.
%    5. Randomly pick exactly XX = round(oxid/100 * N) carbon atoms,  
%         ensuring –OH groups are uniformly distributed.
%    6. For each selected atom:
%         • Assign new C–OH charge/type  
%         • Generate corresponding O and H positions (with randomized orientation)
%         • Update bond and angle connectivity for C–O and O–H links.
%    7. Assemble and inject the new atom/bond/angle blocks into “GrapheneFunct.lt”.
%    8. Clean up all temporary files.
%
% Usage:
%    Adjust the parameters at the top (dir_path, oxid, ...) then
%    run this script. The output “GrapheneFunct.lt” will include the original
%    lattice plus added –OH groups ready for subsequent simulation steps.

%%
clear,clc
tic

%% Settings
dir_path="./";
oxid = 20;                           % Oxidation percentage (e.g 20) 

qCC = 0.265;                            % Charge: atom bound with OH group
qO2 = -0.683;                           % Charge: Oxygen of OH group
qH2 = 0.418;                            % Charge: Hydrogen of OH group

N = 1972;                               % Number of graphene atoms
XX = round(oxid*N/100);              % Number of functionalized atoms

%% Tmp files
graphene_lt = fullfile(dir_path, "graphene.lt");
modified_graphene_lt=fullfile(dir_path, "equilibration/GrapheneFunct.lt");
Bonds_generation = fullfile(dir_path, "DataBonds.txt");
anglesCOH_txt = fullfile(dir_path, "DataAngles.txt");
atomC_txt = fullfile(dir_path, "atomC.txt"); 
atomO_txt = fullfile(dir_path, "atomO.txt"); 
atomH_txt = fullfile(dir_path, "atomH.txt");
atoms_txt = fullfile(dir_path, "atoms.txt");
bondCCO2_txt = fullfile(dir_path, "bondCCO2.txt");
bondO2H2_txt = fullfile(dir_path, "bondO2H2.txt");

%% File reading
if XX==0
    fprintf('Functionalization not required\n')
    fprintf('Creating "%s" as an exact copy of "%s"\n', modified_graphene_lt, graphene_lt)
    return
end
targetWord = '  write("Data Atoms") {';
lineIndex = -1;

try
    copyfile(graphene_lt, modified_graphene_lt);

    fileOrigin = fopen(graphene_lt, 'r');

    lineCounterA = 1;
    while ~feof(fileOrigin)
        line = fgetl(fileOrigin);
        if contains(line, targetWord)
            lineIndex = lineCounterA;
            break;
        end
        lineCounterA = lineCounterA + 1;
    end
    fclose(fileOrigin);
catch Err
    fprintf("%s\n", Err.message)
    return
end


fileOrigin=fopen(graphene_lt,'r');
data=textscan(fileOrigin,'    $atom:CA_%f $mol:m%f @atom:%s %f %f %f %f',...
    'headerlines',lineCounterA);

atomC=data{1,1};
molC=data{1,2};
typeC=string(data{1,3});
chargeC=data{1,4};
xC=data{1,5};
yC=data{1,6};
zC=data{1,7};

%% Bonds generation by distance among atoms

atom_positions = [xC,yC];

k = sqrt(((yC(1)-yC(2)))^2+(xC(1)-xC(2))^2);  % k distance among oblique C atoms
l = abs(yC(3)-yC(2)); % l distance among aligned C atoms
toll=1e-2;

bond_matrix = [];

linestep = size(atom_positions, 1);

for i = 1:linestep
    for j = i+1:linestep
        distance = norm(atom_positions(i, :) - atom_positions(j, :));
        if ( ismembertol(distance, k,toll) || ismembertol(distance, l,toll) )
            bond_matrix = [bond_matrix; [i, j]];
        end
    end
end

%% random atoms functionalization 
attempt=0;
iterations=100; %number of attempts to arrange atoms in a proper way, complying with the distance criterion

% creating a function for the minimum distance between CC atoms (criterion) 
x_known = [20, 99, 197, 394, 592, 789, 968, 1183]; %XX
y_known = [15, 6.49, 4.25, 2.46, 2.456, 1.418, 1.418, 1.418]; %distances

f = @(x) interp1(x_known, y_known, x, 'pchip'); %returns the corresponding minimum spacing between functionalized C atoms


while attempt<=iterations
try
selected_atoms = [];
functPosition = [];
u = 1;
leng=f(XX);
toll = 1e-2;

otheratoms = 1:linestep;

while true

    random_index = randi(length(otheratoms));
    random_atom = otheratoms(random_index);
    x_atom = xC(random_atom);
    y_atom = yC(random_atom);

    if isempty(selected_atoms)
        selected_atoms(u,1) = random_atom;
        functPosition(u, 1) = x_atom;
        functPosition(u, 2) = y_atom;
        u = u + 1;
    else
        valid_distance = true;
        for i = 1:length(selected_atoms)
            for j = 1:length(selected_atoms)
                if i ~= j
                    distance = sqrt((functPosition(i,1) - x_atom)^2 + (functPosition(i,2) - y_atom)^2);
                    if distance <= (leng - toll)
                        valid_distance = false;
                        break;
                    end
                end
            end
            if ~valid_distance
                break;
            end
        end


        if valid_distance 
            selected_atoms(u,1) = random_atom;
            functPosition(u, 1) = x_atom;
            functPosition(u, 2) = y_atom;
            u = u + 1;
        end
    end

    otheratoms(random_index) = [];

    if length(selected_atoms) >= XX
        break; 
    end
end
break
catch 
    attempt=attempt+1;
end

end

if attempt<iterations
for i = 1:XX
    
    j=selected_atoms(i);
    chargeC(j)=qCC;
    typeC(j)="CC";
   
end

fid2=fopen(atomC_txt,'w');
for i=1:linestep
    fprintf(fid2,...
    '    $atom:CA_%1.0f $mol:m%1.0f @atom:%s %1.4f %1.6f %1.6f %1.6f\n',...
    atomC(i),molC(i),typeC(i),chargeC(i),xC(i),yC(i),zC(i));
end
fclose(fid2);

%% Oxygen 

atomO=selected_atoms;
molO=2*ones(XX,1);
typeO=repmat("O2",XX,1);
chargeO=qO2*ones(XX,1);
dist_CC_O2=1.3264; 

for i=1:XX
    index=selected_atoms(i);
    xO(i)=xC(index)';
    yO(i)=yC(index)';
    r=randi(2,1);
    if r==1
        zO(i)=dist_CC_O2; 
    else
        zO(i)=-dist_CC_O2;
    end
end

fid3=fopen(atomO_txt,'w');
for i=1:XX
    fprintf(fid3,...
    '    $atom:O2_%1.0f $mol:m%1.0f @atom:%s %1.4f %1.6f %1.6f %s\n',...
    atomO(i),molO(i),typeO(i),chargeO(i),xO(i),yO(i),zO(i));
end
fclose(fid3);

%% Hydrogen

atomH=selected_atoms;
molH=2*ones(XX,1);
typeH=repmat("H2",XX,1);
chargeH=qH2*ones(XX,1); 
dist_O2_H2=0.945; 
r=dist_O2_H2*cos(23*pi/180); 

for i=1:XX
    index=selected_atoms(i);
    xH(i,1)=xC(index) + r*(-1 + 2*rand()); %23° (113-90)
    
    if rand()<0.5
        yH(i,1)=yC(index) + sqrt(r^2-(xH(i)-xC(index))^2);
    else
        yH(i,1)=yC(index) - sqrt(r^2-(xH(i)-xC(index))^2);
    end

    if zO(i)>0
        zH(i)=zO(i)+sin(dist_O2_H2); 
    else
        zH(i)=zO(i)-sin(dist_O2_H2); 
    end

end

fid4=fopen(atomH_txt,'w');
for i=1:XX
    fprintf(fid4,...
    '    $atom:H2_%1.0f $mol:m%1.0f @atom:%s %1.4f %1.6f %1.6f %s\n',... 
    atomH(i),molH(i),typeH(i),chargeH(i),xH(i),yH(i),zH(i));
end
fclose(fid4);

%% Bonds CC-O2

fid5=fopen(bondCCO2_txt,'w');
for i=1:length(selected_atoms)
    fprintf(fid5,...
    '    $bond:id%1.0f @bond:%s $atom:CA_%1.0f $atom:O2_%1.0f\n',...
    i,"COf",selected_atoms(i),selected_atoms(i));
end
fclose(fid5);

%% Bonds O2-H2
fid6=fopen(bondO2H2_txt,'w');
for i=1:length(selected_atoms)
    fprintf(fid6,...
    '    $bond:id%1.0f @bond:%s $atom:O2_%1.0f $atom:H2_%1.0f\n',...
    i,"OHf",selected_atoms(i),selected_atoms(i));
end
fclose(fid6);

%% Angles

fid7=fopen(anglesCOH_txt,'w');
for i=1:XX
    fprintf(fid7,...
    '    $angle:id%1.0f @angle:COH $atom:CA_%1.0f $atom:O2_%1.0f $atom:H2_%1.0f\n' ,...
    i,selected_atoms(i),selected_atoms(i),selected_atoms(i));
end
fclose(fid7);

%% merging data  

FileContent1 = fileread(atomC_txt);
FileContent2 = fileread(atomO_txt);
FileContent3 = fileread(atomH_txt);
FileContent4 = '  write("Data Bonds") {'; 
FileContent5 = fileread(bondCCO2_txt);
FileContent6 = fileread(bondO2H2_txt);
FileContent7 = '  write("Data Angles") {'; 
FileContent8 = fileread(anglesCOH_txt);

merging1 = [FileContent1, FileContent2, FileContent3];
merging2 = [FileContent5, FileContent6];


fileID = fopen(atoms_txt, 'w');
fprintf(fileID, '%s', merging1);
fprintf(fileID, '%s\n', "}");
fprintf(fileID, '%s\n', FileContent4);
fprintf(fileID, '%s\n', merging2);
fprintf(fileID, '%s\n', "}");
fprintf(fileID, '%s\n', FileContent7);
fprintf(fileID, '%s\n', FileContent8);

fclose(fileID);

%% new functionalized file.lt (GrapheneFunct.lt)

replacementContentA = fileread(atoms_txt);
sourceFileName =modified_graphene_lt;
tempFileName = dir_path+"temp.txt";

replacementStartLineA = lineCounterA+1;
replacementEndLineA = replacementStartLineA + size(atomC, 1)-1;

sourceFile = fopen(sourceFileName, 'r');
tempFile = fopen(tempFileName, 'w');

currentLine = 1;
line = fgets(sourceFile);

while ischar(line)
    if currentLine < replacementStartLineA || currentLine > replacementEndLineA
            fprintf(tempFile, '%s', line);
    elseif currentLine == replacementStartLineA
        fprintf(tempFile, '%s', replacementContentA);
    end
    currentLine = currentLine + 1;
    line = fgets(sourceFile);
end

fclose(sourceFile);
fclose(tempFile);
delete(atoms_txt);
delete(atomC_txt);
delete(atomH_txt);
delete(atomO_txt);
delete(bondCCO2_txt);
delete(bondO2H2_txt);
delete(anglesCOH_txt);
movefile(tempFileName, sourceFileName, 'f');

else 
    fprintf("Number of functionalized atoms not reached\n")
end

fprintf("Elapsed time: %.1f s\n",toc)
