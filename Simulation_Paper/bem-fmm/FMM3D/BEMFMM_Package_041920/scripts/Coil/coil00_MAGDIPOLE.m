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

clear all; %#ok<CLALL>
if ~isunix
    s = pwd; addpath(strcat(s(1:end-5), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-5), '/Engine'));
end

%%  Create the loop structure
%  Define wire cross-section
type = 'circle';
a   = 0;       %   wire radius
M   = 1;       %   number of cross-section subdivisions (individual wires)
%  Define the number of angular subdivisions per loop
N   = 256; 
%   Give the x-coordinates for the centerline of every loop in the xz-plane:
x   = 0.5e-3;         %   centerline intersection with the xz-plane
%   Give the z-coordinates for the centerline of every loop (the z offset)
z  = 0;             %   centerline intersection with the xz-plane
K  = length(z);             %    number of loops in the coil        
%   Create the base model with K loops oriented along the z-axis by cloning one loop with M wires
Pwire1          = [];               % : by 3 point array
Ewire1          = [];               % : by 2 edge array
Swire1          = [];               % : by 1 current strength array for each segment/edge
NFwire1         = ones(K*M, 1);     % : by 1 array of the number of nodes in every wire
S               = N-1;              % number of segments in a single wire 
ringcurrent     = ones(S, 1);       % uniform current weights here
for m = 1:K                 % loop over coil loops   
    [Pwire, Ewire]  = wire_single_loop(x(m), x(m), a, a, N, M, type);  %  this is the major function 
    Pwire(:, 3)     = Pwire(:, 3) + z(m); 
    Ewire1          = [Ewire1; Ewire+size(Pwire1, 1)];
    Pwire1          = [Pwire1; Pwire];
    Swire1          = [Swire1; 1/M*repmat(ringcurrent, M, 1)];
    for n = 1:M
       NFwire1(M*(m-1)+n) = N;  
    end
end

%%   Clone the loop structure
%   Separate the loops along the x-axis
offset = 5e-3;
Pwire2 = Pwire1;
Pwire3 = Pwire1;

Pwire1(:, 1) = Pwire1(:, 1) - offset;
Pwire3(:, 1) = Pwire3(:, 1) + offset;

%   Combine three windings together
Pwire           = [Pwire1; Pwire2; Pwire3];
Ewire           = [Ewire1; Ewire1+size(Pwire1, 1); Ewire1+2*size(Pwire1, 1)];
Swire           = [0.5*Swire1; 0*Swire1; 0.5*Swire1]; %  swap current direction for the second coil (weigfht of -1)
NFwire          = [NFwire1; NFwire1; NFwire1];
segments        = size(Swire, 1); 


%%  Display the coil
Vertices                = size(Pwire, 1);
color(1:Vertices, :)    = repmat([0 0 1], Vertices, 1);
patch('faces', Ewire, 'vertices', Pwire, 'FaceVertexCData', color, 'EdgeColor','interp', 'FaceColor','none','LineWidth', 0.2);
axis 'equal';  axis 'tight';      
xlabel('x, m'); ylabel('y, m'); zlabel('z, m');
view(0, 90);  set(gcf,'Color','White');
title('Coil model (optionally with the relative current strength (-1 to 1))')

%%   Save coil data
strcoil.Pwire   = Pwire;
strcoil.Ewire   = Ewire;
strcoil.Swire   = Swire;
strcoil.NFwire  = NFwire;
save('coil', 'strcoil');
