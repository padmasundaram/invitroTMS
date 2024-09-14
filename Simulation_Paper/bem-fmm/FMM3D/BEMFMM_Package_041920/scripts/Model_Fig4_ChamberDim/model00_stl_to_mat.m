%   This script does stl to MATLAB format conversion. Note that the double
%   points are preserved! Also, performs conversion from mm to m if necessary!

%   Copyright SNM 2017-2019

[~,name,~]=fileparts(pwd); L = length(name)+1;
if ~isunix
    s = pwd; addpath(strcat(s(1:end-L), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-L), '/Engine'));
end

FileName = uigetfile('*.stl','Select the tissue mesh file to open', 'MultiSelect', 'on');
if iscell(FileName)
    for m = 1:length(FileName)
        [P, t, normals, dummy] = stlReadBinary(FileName{m}); 
        P = P*1e-3; %  only if the original data were in mm!
        display(FileName{m});
        NewName =  strcat(FileName{m}(1:end-4), '.mat');
        save(NewName, 'P', 't', 'normals');           
    end
else
        [P, t, normals, dummy] = stlReadBinary(FileName); 
        P = P*1e-3; %  only if the original data were in mm!
        display(FileName);
        NewName =  strcat(FileName(1:end-4), '.mat');
        save(NewName, 'P', 't', 'normals');        
end


