%   This script computes tissue mesh intersections with the YZ plane
%   (saggital plane) located at X=0 or elsewhere, displays and saves the
%   result
%
%   Copyright SNM 2012-2018
%   The Athinoula A. Martinos Center for Biomedical Imaging, Massachusetts General
%   Hospital & ECE Dept., Worcester Polytechnic Inst.

function model03_cross_section_YZ()
    figure();
    [~,name,~]=fileparts(pwd); L = length(name)+1;
    if ~isunix
        s = pwd; addpath(strcat(s(1:end-L), '\Engine'));
    else
        s = pwd; addpath(strcat(s(1:end-L), '/Engine'));
    end

    load CombinedMesh;

    %  Select the YZ-plane window and the x-coordinate of the YZ-plane (all in m)
    ymin = -13e-2;
    ymax = 13e-2;
    zmin = -1e-2;
    zmax = 2e-2;
    X    = 1e-6;  % here is the plane position 

    %   Create coordinates of intersection contours
    tissues = length(name);
    PofYZ = cell(tissues, 1);   %   intersection nodes for a tissue
    EofYZ = cell(tissues, 1);   %   edges formed by intersection nodes for a tissue
    TofYZ = cell(tissues, 1);   %   intersected triangles
    count = [];                 % number of every tissue present in the slice
    for m = 1:tissues
        load(name{m});
        [P, t] = fixmesh(P, t);
        edges = meshconnee(t);
        [Pi, ti, polymask, flag] = meshplaneintYZ(P, t, edges, X);
        if flag % intersection found               
            count               = [count m];
            inslice             = length(count); %  total number of tissues present in the slice
            PofYZ{inslice}      = Pi;           %   intersection nodes
            EofYZ{inslice}      = polymask;     %   edges formed by intersection nodes   
            TofYZ{inslice}      = ti;           %   edges formed by intersection nodes
        end
    end   

    % Display contours
    color = cool(inslice);
    hold on;
    for m = 1:inslice
        for n = 1:size(EofYZ{m}, 1)
            i1 = EofYZ{m}(n, 1);
            i2 = EofYZ{m}(n, 2);
            line(PofYZ{m}([i1 i2], 2), PofYZ{m}([i1 i2], 3), 'Color', color(m, :), 'LineWidth', 2);
        end   
    end

    %   Insert colorbar with tissue names
    mytickmap = cell(inslice, 1);
    for m = 1:inslice
         mytickmap{m} = name{count(m)}(12:end-4);
    end
    ticks      = linspace(0.5/inslice, 1-0.5/inslice, inslice);
    colormap(color);
    colorbar
    colorbar('Ticks', ticks, 'TickLabels', mytickmap);

    title(strcat('Sagittal cross-section at x =', num2str(X), ' m'));
    xlabel('y, m'); ylabel('z, m');
    axis 'equal';  axis 'tight';     
    axis([ymin ymax zmin zmax]);
    set(gcf,'Color','White');

    %   Save the result
    save('YZ_cross_section', 'PofYZ', 'EofYZ', 'TofYZ', 'count', 'color', ...
         'inslice', 'ticks', 'mytickmap', 'ymin', 'ymax', 'zmin', 'zmax', 'X');
end

function [Pi, ti, polymask, flag] = meshplaneintYZ(P, t, edges, X) 
%   This function implements the mesh intersection algorithm for the YZ plane 
%   Output:
%   flag = 0 -> no intersection
%   Pi - intersection nodes
%   polymask - edges formed by intersection nodes

    while 1
        index1  = P(edges(:, 1), 2)==X; 
        index2  = P(edges(:, 2), 2)==X; 
        indexU  = find(index1|index2);
        if ~isempty(indexU)
            X = X*(1+1e-12);
        else
            break;
        end 
    end    
    
    p = [X 0 0];
    n = [1 0 0];
    flag = 1;
    Pi = [];
    ti = [];
    polymask = [];
    
    %   Find all edges (edgesI) intersecting the given plane
    index1  = P(edges(:, 1), 1)>X; 
    index2  = P(edges(:, 2), 1)>X;
    index3  = P(edges(:, 1), 1)<X; 
    index4  = P(edges(:, 2), 1)<X;
    indexU  = index1&index2;
    indexL  = index3&index4;
    indexI  = (~indexU)&(~indexL);
    edgesI  = edges(indexI, :);     %   sorted
    E       = size(edgesI, 1);
    if E==0
        flag = 0;
        return;
    end    
    %   Find all intersection points (Pi)
    N   = repmat(n, [E 1]);                         % repmat plane normal
    Pn  = repmat(p, [E 1]);                         % repmat plane center
    V2      = P(edgesI(:, 2), :);
    V1      = P(edgesI(:, 1), :);
    dot1    = dot(N, (V2 - Pn), 2);
    dot2    = dot(N, (V2 - V1), 2);
    dot3    = dot1./dot2; 
    Pi      = V2 - (V2-V1).*repmat(dot3, [1 3]);
    %  Establish pairs of interconnected edges (pairs of edge nodes)
    AttTriangles   = meshconnet(t, edgesI, 'manifold');    
    AttTrianglesU  = unique([AttTriangles(:, 1)' AttTriangles(:, 2)'])';
    polymask       = zeros(E, 2);
    for m = 1:E
        temp1 = find(AttTriangles(:, 1)==AttTrianglesU(m));
        temp2 = find(AttTriangles(:, 2)==AttTrianglesU(m));
        polymask(m, :)= [temp1 temp2];
    end
    ti = AttTrianglesU;
end
