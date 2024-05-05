%   This is a mesh processor script: it computes basis triangle parameters
%   and necessary potential integrals, and constructs a combined mesh of a
%   multi-object structure (for example, a head or a whole body)
%
%   Copyright SNM/WAW 2017-2020
%clear all %#ok<CLALL>
[~,name,~]=fileparts(pwd); L = length(name)+1;
if ~isunix
    s = pwd; addpath(strcat(s(1:end-L), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-L), '/Engine'));
end

%% Load tissue filenames and tissue display names from index file
index_name = 'tissue_index.txt';
[name, tissue, cond, enclosingTissueIdx] = tissue_index_read(index_name);

%%  Generic tissue list (for graphics only)
% tissue{1} = ' scalp'; 
% tissue{2} = ' skull';
% tissue{3} = ' CSF'; 
% tissue{4} = ' GM';
% tissue{5} = ' cerebellum'; 
% tissue{6} = ' WM'; 
% tissue{7} = ' ventricles';

%%  Load tissue meshes and combine individual meshes into a single mesh
tic
PP = [];
tt = [];
nnormals = [];
Indicator = [];

%   Combine individual meshes into a single mesh
for m = 1:length(name)
    load(name{m}); 
    tt = [tt; t+size(PP, 1)];
    PP = [PP; P];
    nnormals = [nnormals; normals];    
    Indicator= [Indicator; repmat(m, size(t, 1), 1)];
    disp(['Successfully loaded file [' name{m} ']']);
end
t = tt;
P = PP;
normals = nnormals;
LoadBaseDataTime = toc

%%  Fix triangle orientation (just in case, optional)
tic
t = meshreorient(P, t, normals);
%%   Process other mesh data
Center      = 1/3*(P(t(:, 1), :) + P(t(:, 2), :) + P(t(:, 3), :));  %   face centers
Area        = meshareas(P, t);  
SurfaceNormalTime = toc

%%  Assign facet conductivity information
tic
condambient = 0.0; %   air
[contrast, condin, condout] = assign_initial_conductivities(cond, condambient, Indicator, enclosingTissueIdx);
InitialConductivityAssignmentTime = toc

%%  Check for and process triangles that have coincident centroids
tic
disp('Checking combined mesh for duplicate facets ...');
[P, t, normals, Center, Area, Indicator, condin, condout, contrast] = ...
    clean_coincident_facets(P, t, normals, Center, Area, Indicator, condin, condout, contrast);
disp('Resolved all duplicate facets');
N           = size(t, 1);
DuplicateFacetTime = toc

%%   Find topological neighbors
tic
DT = triangulation(t, P); 
tneighbor = neighbors(DT);
% Fix cases where not all triangles have three neighbors
tneighbor = pad_neighbor_triangles(tneighbor);

%%   Save base data
FileName = 'CombinedMesh.mat';
save(FileName, 'P', 't', 'normals', 'Area', 'Center', 'Indicator', 'name', 'tissue', 'cond', 'enclosingTissueIdx', 'condin', 'condout', 'contrast');
ProcessBaseDataTime = toc

%%   Add accurate integration for electric field/electric potential on neighbor facets
%   Indexes into neighbor triangles
numThreads = 18;        %   number of cores to be used
RnumberE        = 4;    %   number of neighbor triangles for analytical integration (fixed, optimized)
ineighborE      = knnsearch(Center, Center, 'k', RnumberE);   % [1:N, 1:Rnumber]
ineighborE      = ineighborE';           %   do transpose    
EC         = meshneighborints_2(P, t, normals, Area, Center, RnumberE, ineighborE, numThreads);

%%   Normalize sparse matrix EC by variable contrast (for speed up)
N   = size(Center, 1);
ii  = ineighborE;
jj  = repmat(1:N, RnumberE, 1); 
CO  = sparse(ii, jj, contrast(ineighborE));
EC  = CO.*EC;

tic
NewName  = 'CombinedMeshP.mat';
save(NewName, 'tneighbor',  'RnumberE',   'ineighborE', 'EC', '-v7.3');
SaveBigDataTime = toc



%%
% %   This is a mesh processor script: it computes basis triangle parameters
% %   and necessary potential integrals, and constructs a combined mesh of a
% %   multi-tissue structure (a head or a whole body)
% %
% %   Copyright SNM 2012-2018
% %   The Athinoula A. Martinos Center for Biomedical Imaging, Massachusetts General
% %   Hospital & ECE Dept., Worcester Polytechnic Inst.
% 
% % clear all %#ok<CLALL>
% if ~isunix
%     s = pwd; addpath(strcat(s(1:end-5), '\Engine'));
% else
%     s = pwd; addpath(strcat(s(1:end-5), '/Engine'));
% end
% 
% Rnumber   = 12;     %   number of neighbor triangles for analytical integration (fixed, optimized)  
% %    Replaces R: dimensionless radius of an enclosing sphere for precise 
% %                integration
% %    Rnumber=12 - integrals for 12 closest neighbor triangles are calculated
% %    precisely and stored, etc.
% 
% gauss = 25; %    Number of integration points in the Gaussian quadrature 
% %                for the outer potential integrals
%             %    Numbers 1, 4, 7, 13, 25 are permitted 
% eps0  = 8.85418782e-012;  %   Dielectric permittivity of vacuum(~air)   
% 
% %%  Generic tissue list (for graphics)
% tissue{1} = ' cerebellum'; 
% tissue{2} = ' bath';
% tissue{3} = ' chamber';
% 
% %%  Load tissue meshes and combine individual meshes into a single mesh
% tic
% PP = [];
% tt = [];
% nnormals = [];
% Indicator = [];
% 
% %   Compile tissue list (cell aray)
% name{1}       = 'head_brain.mat';            %     file to import
% name{2}       = 'head_chin.mat';           %     file to import
% name{3}       = 'head_chout.mat';           %     file to import
% 
% 
% 
% 
% % Remove somew tissues if necessary
% remove          = [];
% name(remove)    = [];
% 
% %   Combine individual meshes into a single mesh
% for m = 1:length(name)
%     load(name{m});    
%     tt = [tt; t+size(PP, 1)];
%     PP = [PP; P]; % in m !
%     nnormals = [nnormals; normals];    
%     Indicator= [Indicator; repmat(m, size(t, 1), 1)];
% end
% t = tt;
% P = PP;
% normals = nnormals;
% LoadBaseDataTime = toc
% 
% %%  Process base mesh data
% %  Process surface normal data
% tic
% N           = size(t, 1);
% for m = 1:N
%     Vertexes        = P(t(m, 1:3)', :)';
%     r1              = Vertexes(:, 1);
%     r2              = Vertexes(:, 2);
%     r3              = Vertexes(:, 3);
%     tempv           = cross(r2-r1, r3-r1);  %   definition (*)
%     temps           = sqrt(tempv(1)^2 + tempv(2)^2 + tempv(3)^2);
%     normalcheck     = tempv'/temps;
%     if sum(normalcheck.*normals(m, :))<0;   %   rearrange vertices to have exactly the outer normal
%         t(m, 2:3) = t(m, 3:-1:2);           %   by definition (*)
%     end     
% end   
% %   Process other data
% Center      = zeros(N, 3);   %   face center
% Area        = zeros(N, 1);   %   face area
% Size        = zeros(N, 1);   %   face size       
% ineighbor   = cell(N, 1);    %   array of neighbor triangles
% Center      = 1/3*(P(t(:, 1), :) + P(t(:, 2), :) + P(t(:, 3), :)); 
% Area        = meshareas(P, t);
% Size        = sqrt(mean(Area));
% %   Save base data
% FileName = 'CombinedMesh.mat';
% save(FileName, 'P', 't', 'normals', 'Area', 'Center', 'Indicator', 'Size', 'Rnumber', 'name', 'remove', 'tissue');
% ProcessBaseDataTime = toc
% 
% %%   Add accurate integration 
% tic
% %   Indexes into neighbor triangles
% ineighbor  = knnsearch(Center, Center, 'k', Rnumber);   % [1:N, 1:Rnumber]
% 
% iintegralx   = zeros(Rnumber, N);    %   integral component for array of neighbor triangles
% iintegraly   = zeros(Rnumber, N);    %   integral component for array of neighbor triangles
% iintegralz   = zeros(Rnumber, N);    %   integral component for array of neighbor triangles
% DoTriangleNeighborSeearchTime = toc
% 
% %   Gaussian weights for analytical integration (for the outer integral)
% if gauss == 1;  [coeffS, weightsS, IndexS]  = tri(1, 1); end;
% if gauss == 4;  [coeffS, weightsS, IndexS]  = tri(4, 3); end;
% if gauss == 7;  [coeffS, weightsS, IndexS]  = tri(7, 5); end;
% if gauss == 13; [coeffS, weightsS, IndexS]  = tri(13, 7); end;
% if gauss == 25; [coeffS, weightsS, IndexS]  = tri(25, 10); end;
% W           = repmat(weightsS', 1, 3);
% 
% %   Main loop for analytical double integrals (parallel, 24 workers)
% tic
% %parpool(24);
% parfor m = 1:N     
%     r1      = P(t(m, 1), :);    %   row
%     r2      = P(t(m, 2), :);    %   row
%     r3      = P(t(m, 3), :);    %   row   
%     index   = ineighbor(m, :);      
%     ObsPoints   = zeros(Rnumber*IndexS, 3);  
%     I           = zeros(Rnumber, 3);
%     for n = 1:Rnumber
%         num = index(n);
%         for p = 1:IndexS
%             ObsPoints(p+(n-1)*IndexS, :)  = coeffS(1, p)*P(t(num, 1), :) +  coeffS(2, p)*P(t(num, 2), :) +  coeffS(3, p)*P(t(num, 3), :);
%         end
%     end    
%     J = potint2(r1, r2, r3, normals(m, :), ObsPoints); %   I was calculated with the area Area(n)  
%     for n = 1:Rnumber
%         I(n, :) = sum(W.*J([1:IndexS]+(n-1)*IndexS, :), 1); 
%     end
%     iintegralx(:, m) = I(:, 1);    %   accurate integrals
%     iintegraly(:, m) = I(:, 2);    %   accurate integrals
%     iintegralz(:, m) = I(:, 3);    %   accurate integrals     
%     %   Subtract the central-point contribution in advance and redifine integrals
%     const = 1/(4*pi*eps0);                         
%     k     = index==m;
%     temp    = repmat(Center(m, :), Rnumber, 1) - Center(index, :); %   these are distances to the observation/target triangle
%     DIST    = sqrt(dot(temp, temp, 2));                            %   single column                
%     I       = Area(m)*temp./repmat(DIST.^3, 1, 3);                 %   center-point integral, standard format                                                                                                 
%     contributionx = const*(I(:, 1) - iintegralx(:, m));
%     contributiony = const*(I(:, 2) - iintegraly(:, m));
%     contributionz = const*(I(:, 3) - iintegralz(:, m));
%     contributionx(k) = 0;
%     contributiony(k) = 0;
%     contributionz(k) = 0;          
%     iintegralx(:, m) = contributionx;
%     iintegraly(:, m) = contributiony;
%     iintegralz(:, m) = contributionz;    
% end 
% MainLoopParallelTime = toc
% 
% tic
% NewName  = 'CombinedMeshP.mat';
% save(NewName, 'Rnumber', 'Size', 'gauss', 'ineighbor', 'iintegralx', 'iintegraly', 'iintegralz');
% SaveBigDataTime = toc
% delete(gcp('nocreate'));
