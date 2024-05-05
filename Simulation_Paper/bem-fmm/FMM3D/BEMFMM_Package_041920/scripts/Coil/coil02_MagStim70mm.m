%   "Wire coil creator" tool for the entire coil. Models the outer surface
%   of coil windings by multiple wire segments. The output is saved in the
%   binary file coil.mat and includes:
%   Pwire{:, 3} - set of nodes for all wires 
%   Ewire{:, 2} - set of edges/segments for all wires (current flows from the first
%   edge node to the second edge node)
%   Swire{:, 1} - current strength weight for every segment/edge (from -1 to 1)
%   This parameter depends on the number of wires M in the indivudual
%   conductor and is scaled by default as 1/M
%   NFwire(:, 1): number of nodes in every individual wire (optional)
%   PFwire(:, 1): set of nodes in every individual wire (optional)
%   This is a cell array cell(:, 1) used for graphics
%
%   Copyright SNM 2017-2018
%   The Athinoula A. Martinos Center for Biomedical Imaging, Massachusetts General
%   Hospital & ECE Dept., Worcester Polytechnic Inst.

%---------------------------------------------------------
% MagStim 70 mm Figure of 8 coil 
% Inside Diameter 56 mm
% Outside Diameter 87 mm
% number of turns 9
% Inductance 15.5 uH
% 2.2 T peak magnetic field
% 660 V/m peak E-field
% Flat spiral winding
% WC = 1.75 mm x 6 mm (wire cross section)
%---------------------------------------------------------


clear all; %#ok<CLALL>
if ~isunix
    s = pwd; addpath(strcat(s(1:end-5), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-5), '/Engine'));
end

%%  Create the loop structure
%  Define wire cross-section
type = 'rect';
a   = 1.75e-3;       %   wire width
b = 6e-3; % wire length
M   = 6; %;                   %   number of cross-section subdivisions (individual wires)
%  Define the number of angular subdivisions per loop
N   = 64; 

flag = 2;
sk   = 0;           %   uniform current distribution (skin layer)



%   Create the base model with K loops oriented along the z-axis by cloning one loop with M wires
strcoil.Swire       = [];
strcoil.Ewire       = [];
strcoil.Pwire       = [];
S               = N-1;              % number of segments in a single wire 

P                   = [];
t                   = [];
tind                = [];


%   Creating the coil turns in a loop
offset = 45e-3;
a0 = [28.8750,   30.5938,   32.3125,   34.0312,   35.7500,   37.4688,   39.1875,   40.9062,   42.6250]*1e-3;
K  = length(a0);             %    number of loops in the coil      
for m = 1:K                 % loop over coil loops   
   
    theta   = linspace(0, 2*pi, 100);
    
    y = a0(m)*cos(theta);                 %   planar ellipse
    x = a0(m)*sin(theta);                 %   planar ellipse
    z = zeros(size(x));
    
    %[Pwire, Ewire]  = wire_single_loop(y(m), x(m), a, b, N, M, type);  %  this is the major function 
    %   Create CAD and wire models for the single conductor
    Pcenter(:, 1) = x'+ offset;
    Pcenter(:, 2) = y' ;
    Pcenter(:, 3) = z';
    
    strcoil_temp    = meshwire(Pcenter, a, b, M, flag, sk);
    [Ptemp, ttemp]  = meshsurface(Pcenter, a, b, M, flag);
    [Ptemp, ttemp]  = meshfix(Ptemp, ttemp);  
    tindtemp       = 1*ones(size(t, 1), 1);
    
    %   Accumulate loops 
    strcoil.Swire       = [strcoil.Swire; strcoil_temp.Swire];
    strcoil.Ewire       = [strcoil.Ewire; strcoil_temp.Ewire+size(strcoil.Pwire, 1)];
    strcoil.Pwire       = [strcoil.Pwire; strcoil_temp.Pwire];     

    t          = [t; ttemp+size(P, 1)];    
    tind       = [tind; tindtemp];    
    P          = [P; Ptemp];  
end

%   Creating the coil turns in a loop
offset = -45e-3;
a0 = -[28.8750,   30.5938,   32.3125,   34.0312,   35.7500,   37.4688,   39.1875,   40.9062,   42.6250]*1e-3;
K  = length(a0);             %    number of loops in the coil      
for m = 1:K                 % loop over coil loops   
   
    theta   = linspace(0, 2*pi, 100);
    
    y = a0(m)*cos(theta);                 %   planar ellipse
    x = a0(m)*sin(theta);                 %   planar ellipse
    z = zeros(size(x));
    
    %[Pwire, Ewire]  = wire_single_loop(y(m), x(m), a, b, N, M, type);  %  this is the major function 
    %   Create CAD and wire models for the single conductor
    Pcenter(:, 1) = x'+ offset;
    Pcenter(:, 2) = y' ;
    Pcenter(:, 3) = z';
    
    strcoil_temp    = meshwire(Pcenter, a, b, M, flag, sk);
    [Ptemp, ttemp]  = meshsurface(Pcenter, a, b, M, flag);
    [Ptemp, ttemp]  = meshfix(Ptemp, ttemp);  
    tindtemp       = 1*ones(size(t, 1), 1);
    
    %   Accumulate loops 
    strcoil.Swire       = [strcoil.Swire; -strcoil_temp.Swire];
    strcoil.Ewire       = [strcoil.Ewire; strcoil_temp.Ewire+size(strcoil.Pwire, 1)];
    strcoil.Pwire       = [strcoil.Pwire; strcoil_temp.Pwire];     

    t          = [t; ttemp+size(P, 1)];    
    tind       = [tind; tindtemp];    
    P          = [P; Ptemp];  
end



%% display
%   Display CAD and wire models
f1 = figure;
bemf1_graphics_coil_CAD(P, t, 0);
view(0,30);


 
 
%%   Save coil data
save('coil', 'strcoil');
save('coilCAD', 'P', 't', 'tind');  %   optional, slow