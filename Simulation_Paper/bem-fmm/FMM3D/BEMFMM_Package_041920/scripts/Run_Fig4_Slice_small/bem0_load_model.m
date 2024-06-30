%   This script loads mesh data into the MATLAB workspace. The data include
%   surface meshes and the potential integrals. It also computes useful
%   sparse matrices for GMRES speed up
%
%   Copyright SNM/WAW 2017-2020

%%  Add path (Windows/Linux commands are different)

s = pwd;
if(~isunix)
    slash = '\..\';
else
    slash = '/../';
end

%   If a simulation has already been run with a different model,
%   that model's path may have been loaded, and the desired model's 
%   files would then be shadowed by the files on the old model path
warning off; rmpath(genpath(s)); warning on;

engine_path =   [s, slash, 'Engine']; 
model_path =    [s, slash, 'Model_Fig3_Slice_small']; %  <------- ***** This line controls which model is loaded
coil_path =     [s, slash, 'Coil'];

addpath(engine_path);
addpath(model_path);
addpath(coil_path);

%%  Define EM constants
eps0        = 8.85418782e-012;  %   Dielectric permittivity of vacuum(~air)
mu0         = 1.25663706e-006;  %   Magnetic permeability of vacuum(~air)

%%  Import geometry and electrode data. Create useful sparse matrices. Import existing solution (if any)
tic
h                   = waitbar(0.5, 'Please wait - loading model data and existing solution (if any)');
load CombinedMesh;
load CombinedMeshP;


close(h);
LoadTime = toc
