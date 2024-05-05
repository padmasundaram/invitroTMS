%   This script does pcode conversion. 

%   Copyright SNM 2017-2019

clear all; %#ok<CLALL>

FileName = uigetfile('*.m','Select the tissue mesh file to open', 'MultiSelect', 'on');
if iscell(FileName)
    for m = 1:length(FileName)        
        display(FileName{m});
        pcode(FileName{m});
    end
else        
        display(FileName);
        pcode(FileName);
end