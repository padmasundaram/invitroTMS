%   This script creates the mesh (both CAD surface mesh and a computational
%   wire grid) for a MagVenture figure eight planar coil with 1 A of total
%   current
%   The output is saved in the binary file coil.mat and includes:
%   strcoil.Pwire(:, 3) - set of nodes for all wires 
%   strcoil.Ewire(:, 2) - set of edges or current dipoles for all wires
%   (current flows from the first edge node to the second edge node)
%   strcoil.Swire{:, 1} - current strength weight for every elementary
%   dipole asssuring that the total conductor current through any
%   cross-section is 1 A.
%
%   Copyright SNM 2018-2019

clear all; 
if ~isunix
    s = pwd; addpath(strcat(s(1:end-5), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-5), '/Engine'));
end

%   The coil is in the form of two interconnected spiral arms. The
%   conductor centerline model is given first
theta = [3*pi/2:pi/25:12*pi];
a0 = 0.017; b0 = 0.0006;
r = a0 + b0*theta;                  %   Archimedean spiral
x1 = r.*cos(theta);                 %   first half
y1 = r.*sin(theta);                 %   first half
x2 = 2*x1(end) - x1(end-1:-1:1);    %   second half
y2 = 2*y1(end) - y1(end-1:-1:1);    %   second half
x = [x1 x2];  y = [y1 y2];          %   join both halves
x = x - mean(x);                    %   center the curve

% plot(x, y, '*-'); axis equal; grid on; title('Conductor centerline')

%   Other parameters
a    = 3.60e-3;     %   z-side, m  (for a rectangle cross-section)
b    = 2.61e-3;     %   x-side, m  (for a rectangle cross-section)
M    = 20;          %   number of cross-section subdivisions 
flag = 2;           %   rect. cross-section    
sk   = 1;           %   surface current distribution (skin layer)

%   Create CAD and wire models for the single conductor
Pcenter(:, 1) = x';
Pcenter(:, 2) = y';
Pcenter(:, 3) = a/2;
strcoil       = meshwire(Pcenter, a, b, M, flag, sk);
[P, t]        = meshsurface(Pcenter, a, b, M, flag);  %   CAD mesh (optional, slow)     
[P, t]        = meshfix(P, t);  
tind          = 1*ones(size(t, 1), 1);

Ewire       = [];
Pwire       = []; 
Swire       = [];
Pa          = [];
ta          = []; 

%   Construct two CAD and wire models 
strcoil.Swire       = [strcoil.Swire; strcoil.Swire];
strcoil.Ewire       = [strcoil.Ewire; strcoil.Ewire+size(strcoil.Pwire, 1)];
Pwire               = strcoil.Pwire;
Pwire(:, 3)         = Pwire(:, 3) + 9.2e-3;
strcoil.Pwire       = [strcoil.Pwire; Pwire]; 
strcoil.Pwire(:, 3) = strcoil.Pwire(:, 3) - min(strcoil.Pwire(:, 3));

t          = [t; t+size(P, 1)];
tind       = [tind; 2*tind];
Pup        = P;
Pup(:, 3)  = Pup(:, 3) + 9.2e-3;
P          = [P; Pup];
P(:, 3)    = P(:, 3) - min(P(:, 3)); 

%   Display CAD and wire models
bemf1_graphics_coil_CAD(P, t, 0);
%bemf1_graphics_coil_wire(strcoil, [0 1 0]); 
view(0, 0);
 
save('coil', 'strcoil');
save('coilCAD', 'P', 't', 'tind');  %   optional, slow
