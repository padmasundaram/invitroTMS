%   This script assigns coil position, dIdt (for
%   electric field), I0 (for magnetic field), and displays the resulting
%   head-coil geometry. It may be be executed several times to adjust the
%   coil position
%
%   Copyright SNM/WAW 2017-2020

addpath('../../../../FMM3D/matlab');

%%  Define dIdt (for electric field)

%%  Coil parameters
%  Define dIdt (for electric field)
dIdt = 5e6;       %   Amperes/sec (2*pi*I0/period)
pw = 4.3e-6; % initial ramp time

%  Define I0 (for magnetic field)
I0 = pw * dIdt; %5e3;       %   Amperes

%%  Define field margin (for plotting)
margin      = 0.80;     %   Only for fields plotting

%%  Load base coil data, define coil excitation/position, define coil array if necesary
if exist('strcoil', 'var')
    clear strcoil;    
end
load coil.mat;       

%%   Define coil position: rotate and then tilt and move the entire coil as appropriate
coilaxis        = [0 0 1];                  %   Transformation 1: rotation axis
theta           = 0*pi/180.0;           %   Transformation 1: angle to rotate about axis
Nx = 0; Ny = +0.0; Nz = -1;             %   Transformation 2: New coil centerline direction

%MoveX = 0; MoveY = -12e-2; MoveZ = 5.5e-2; %   Transformation 3: New coil position
MoveX = 0; MoveY = 0.0; MoveZ = 0.5e-3 -3e-3; %   Transformation 3: New coil position

%   Apply Transformation 1: rotation about coil centerline
strcoil.Pwire   = meshrotate2(strcoil.Pwire, coilaxis, theta);
% Coil.P          = meshrotate2(Coil.P, coilaxis, theta);

%   Apply Transformation 2: Tilt the coil axis with direction vector Nx, Ny, Nz as required
strcoil.Pwire = meshrotate1(strcoil.Pwire, Nx, Ny, Nz);
% Coil.P        = meshrotate1(Coil.P, Nx, Ny, Nz);

%   Apply Transformation 3: Move the coil as required
strcoil.Pwire(:, 1)   = strcoil.Pwire(:, 1) + MoveX;
strcoil.Pwire(:, 2)   = strcoil.Pwire(:, 2) + MoveY;
strcoil.Pwire(:, 3)   = strcoil.Pwire(:, 3) + MoveZ;

% Coil.P(:, 1)   = Coil.P(:, 1) + MoveX;
% Coil.P(:, 2)   = Coil.P(:, 2) + MoveY;
% Coil.P(:, 3)   = Coil.P(:, 3) + MoveZ;

%%   Define the observation line from the bottom center of the coil into the head
M = 10000;        
argline      = linspace(0, 10e-3, M);              %   distance along a 100 mm long line   
dirline      = -[Nx Ny Nz]/norm([Nx Ny Nz]);        %   line direction (along the coil axis)   
offline      = 0e-3;                                %   offset from the coil
pointsline(1:M, 1) = MoveX + dirline(1)*(argline + offline);
pointsline(1:M, 2) = MoveY + dirline(2)*(argline + offline);
pointsline(1:M, 3) = MoveZ + dirline(3)*(argline + offline);

save('output_coil_data', 'strcoil', 'argline', 'pointsline');

%%  Head graphics
figure('color','w');
tissue_to_plot = 'Ringer';
t0 = t(Indicator==find(strcmp(tissue, tissue_to_plot)), :);    % (change indicator if necessary: 1-skin, 2-skull, etc.)
str.EdgeColor = 'k'; str.FaceColor = [0.5 0.5 0.5]; str.FaceAlpha = 0.3; 
bemf2_graphics_base(P, t0, str);

%% Coil graphics    
hold on;
% bemf1_graphics_coil_CAD(Coil.P, Coil.t, 0);
bemf1_graphics_coil_wire(strcoil, 'r');

%%  Plot coil centerline
hold on;
plot3(pointsline(:, 1), pointsline(:, 2), pointsline(:, 3), '-r', 'lineWidth', 3);

%%  General settings
axis 'equal';  axis 'tight';   
daspect([1 1 1]);
set(gcf,'Color','White');
lighting phong;
view(180, 90); axis off; camzoom(1.0);
set(gcf,"Position",[81 431 1014 1017]);

%% Find nearest intersections of the coil centerline w tissues
%   Ray parameters (in mm here)
orig = 1e3*[MoveX MoveY MoveZ];     %   ray origin
dir  = dirline;     %   ray direction
dist = 10000;        %   ray length (finite segment, in mm here)

intersections_to_find = tissue;

for m = 1:length(intersections_to_find)
    k = find(strcmp(intersections_to_find{m}, tissue));
    disp(intersections_to_find{m});
    S = load(name{k});
    
    d = meshsegtrintersection(orig, dir, dist, S.P, S.t);
    IntersectionPoint = min(d(d>0))
    if ~isempty(IntersectionPoint)
        Position = orig + dir*IntersectionPoint
    end
    sprintf(newline);
end
