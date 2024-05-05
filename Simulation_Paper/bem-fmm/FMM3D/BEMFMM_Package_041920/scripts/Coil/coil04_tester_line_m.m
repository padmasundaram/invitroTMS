%   "Coil tester line magnetic" tool for the entire coil. Outputs
%   total magnetic field (or any of the field components) along the line given I0 A
%   of conductor current

%   Copyright SNM 2017-2019

clear all; %#ok<CLALL>
if ~isunix
    s = pwd; addpath(strcat(s(1:end-5), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-5), '/Engine'));
end
load coil; load coilCAD;
addpath('/cluster/fusion/padma/bem-fmm/teppei_fmm3d_092022/FMM3D/matlab/');

%%  Define EM constants
mu0         = 1.25663706e-006;  %   Magnetic permeability of vacuum(~air)

%%  Coil parameters
%  Define dIdt (for electric field)
dIdt = 9.4e7;       %   Amperes/sec (2*pi*I0/period)
%  Define I0 (for magnetic field)
I0 = 5e3;       %   Amperes

%%  Define nodal points along the line (1xM nodal points)
M = 1001;        
line      = linspace(-0.05, 0.05, M);
pointsline(1:M, 1) = 0;
pointsline(1:M, 2) = 0;
pointsline(1:M, 3) = line';   

%%   Plot the line
f1 = figure;
hold on;
bemf1_graphics_coil_CAD(P, t, 0);
plot3(pointsline(:, 1), pointsline(:, 2), pointsline(:, 3), '-r', 'lineWidth', 2);
view(10, 20);

%%  Find the B-field on the line 
tic
Binc          = bemf3_inc_field_magnetic(strcoil, pointsline, I0, mu0); 
fieldTime = toc

%%  Graphics for the line
f2 = figure;
hold on;
plot(line, +I0*Binc(:, 1), '-r', 'LineWidth', 2); 
plot(line, +I0*Binc(:, 2), '-m', 'LineWidth', 2); 
plot(line, +I0*Binc(:, 3), '-b', 'LineWidth', 2); 
grid on;

title('Line field B in T, red:Bx magenta:By blue:Bz');
xlabel('Distance z from the origin, mm');
ylabel('Magnetic field, T');
set(gcf,'Color','White');