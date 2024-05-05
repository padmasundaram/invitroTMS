


%%  Define all three planes (m)

X       = 1e-6;                 %   YZ Cross-section position, m
Y       = 1e-6;                     %   XZ Cross-section position, m
Z       = 200e-6;                 %   XY Cross-section position, m
delta   = 1.5e-3;                       %   half plane window, m
xmin = - delta;                   % Cross-section left edge
xmax = + delta;                   % Cross-section right edge
ymin = - delta;                   % Cross-section posterior edge
ymax = + delta;                   % Cross-section anterior edge
yminzoom = -1.5e-3;                   % Cross-section posterior edge
ymaxzoom = +1.5e-3;                   % Cross-section anterior edge
zmin = -0.25e-3;                % Cross-section inferior edge
zmax = 0.4e-3 + 0.25e-3;        % Cross-section superior edge

%%   Process cross-section data to enable fast (real time) display 
%   This block finds all edges and attached triangles for separate brain
%   compartments. This script is required for subsequent visualizations.
%   Process surface model data
tic
%   Preallocate cell arrays
m_max = length(tissue);
tS = cell(m_max, 1);
nS = tS; %  Reuse this empty cell array for other initialization
eS = tS;
TriPS = tS;
TriMS = tS;
ENinside = tS;
ENoutside = tS;
PS = P; % * 1e3; % Convert to mm
for m = 1:m_max
    tS{m} = t(Indicator == m, :);
    nS{m} = normals(Indicator == m, :);
    [eS{m}, TriPS{m}, TriMS{m}] = mt(tS{m}); 
    ENinside{m}  = Eninside(Indicator == m);
    ENoutside{m} = Enoutside(Indicator == m);
end
SurfaceDataProcessTime = toc

%%  Figure with planes
figure;
%%  Plot the planes
patch([xmin xmin xmax xmax],[ymin ymax ymax ymin], [Z Z Z Z], 'c', 'FaceAlpha', 0.35);
patch([xmin xmin xmax xmax],[Y Y Y Y], [zmin zmax zmax zmin], 'c', 'FaceAlpha', 0.35);
patch([X X X X], [ymin ymin ymax ymax], [zmin zmax zmax zmin], 'c', 'FaceAlpha', 0.35);

%%  Head graphics
tissue_to_plot = 'Ringer';
t0 = t(Indicator==find(strcmp(tissue, tissue_to_plot)), :);    % (change indicator if necessary: 1-skin, 2-skull, etc.)
str.EdgeColor = 'none'; str.FaceColor = [1 0.75 0.65]; str.FaceAlpha = 0.2; 
bemf2_graphics_base(P, t0, str);
tissue_to_plot = 'Slice';
t0 = t(Indicator==find(strcmp(tissue, tissue_to_plot)), :);    % (change indicator if necessary: 1-skin, 2-skull, etc.)
str.EdgeColor = 'none'; str.FaceColor = [1 0.75 0.65]; str.FaceAlpha = 1; 
bemf2_graphics_base(P, t0, str);
%% Coil graphics    
hold on;
bemf1_graphics_coil_wire(strcoil, 'k');

%%  General settings
axis 'equal';  axis 'tight';   
daspect([1 1 1]);
set(gcf,'Color','White');
camlight; lighting phong;
view(157, 25); camzoom(2)

%%  XY plane
figure;

tissues = length(name);
PofXY = cell(tissues, 1);   %   intersection nodes for a tissue
EofXY = cell(tissues, 1);   %   edges formed by intersection nodes for a tissue
TofXY = cell(tissues, 1);   %   intersected triangles
NofXY = cell(tissues, 1);   %   normal vectors of intersected triangles
countXY = [];   %   number of every tissue present in the slice
for m = 1:tissues 
    [Pi, ti, polymask, flag] = meshplaneintXY(PS, tS{m}, eS{m}, TriPS{m}, TriMS{m}, Z);
    if flag % intersection found                
        countXY               = [countXY m];
        PofXY{m}            = Pi;               %   intersection nodes
        EofXY{m}            = polymask;         %   edges formed by intersection nodes
        TofXY{m}            = ti;               %   intersected triangles
        NofXY{m}            = nS{m}(ti, :);     %   normal vectors of intersected triangles        
    end
end
color   = prism(length(tissue)); color(4, :) = [0 1 1];
for m = countXY
    edges           = EofXY{m};             %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXY{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXY{m}(:, 2);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
patch([xmin xmin xmax xmax],[ymin ymax ymax ymin], 'c', 'FaceAlpha', 0.35);
xlabel('x, m'); ylabel('y, m');
axis 'equal';  axis 'tight'; 
set(gcf,'Color','White');    

%%  XZ plane
figure;

tissues = length(name);
PofXZ = cell(tissues, 1);   %   intersection nodes for a tissue
EofXZ = cell(tissues, 1);   %   edges formed by intersection nodes for a tissue
TofXZ = cell(tissues, 1);   %   intersected triangles
NofXZ = cell(tissues, 1);   %   normal vectors of intersected triangles
countXZ = [];   %   number of every tissue present in the slice
for m = 1:tissues 
    [Pi, ti, polymask, flag] = meshplaneintXZ(PS, tS{m}, eS{m}, TriPS{m}, TriMS{m}, Y);
    if flag % intersection found                
        countXZ               = [countXZ m];
        PofXZ{m}            = Pi;               %   intersection nodes
        EofXZ{m}            = polymask;         %   edges formed by intersection nodes
        TofXZ{m}            = ti;               %   intersected triangles
        NofXZ{m}            = nS{m}(ti, :);     %   normal vectors of intersected triangles        
    end
end
color   = prism(length(tissue)); color(4, :) = [0 1 1];
for m = countXZ
    edges           = EofXZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXZ{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
patch([xmin xmin xmax xmax], [zmin zmax zmax zmin], 'c', 'FaceAlpha', 0.35);
xlabel('x, m'); ylabel('z, m');
axis 'equal';  axis 'tight'; 
set(gcf,'Color','White');

%%  YZ plane

figure;

tissues = length(name);
PofYZ = cell(tissues, 1);   %   intersection nodes for a tissue
EofYZ = cell(tissues, 1);   %   edges formed by intersection nodes for a tissue
TofYZ = cell(tissues, 1);   %   intersected triangles
NofYZ = cell(tissues, 1);   %   normal vectors of intersected triangles
countYZ = [];   %   number of every tissue present in the slice
for m = 1:tissues 
    [Pi, ti, polymask, flag] = meshplaneintYZ(PS, tS{m}, eS{m}, TriPS{m}, TriMS{m}, X);
    if flag % intersection found                
        countYZ               = [countYZ m];
        PofYZ{m}            = Pi;               %   intersection nodes
        EofYZ{m}            = polymask;         %   edges formed by intersection nodes
        TofYZ{m}            = ti;               %   intersected triangles
        NofYZ{m}            = nS{m}(ti, :);     %   normal vectors of intersected triangles        
    end
end
color   = prism(length(tissue)); color(4, :) = [0 1 1];
for m = countYZ
    edges           = EofYZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofYZ{m}(:, 2);       %   this is for the contour  
    points(:, 2)    = +PofYZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
patch([ymin ymin ymax ymax], [zmin zmax zmax zmin], 'c', 'FaceAlpha', 0.35);
xlabel('y, m'); ylabel('z, m');
axis 'equal';  axis 'tight'; 
set(gcf,'Color','White');

