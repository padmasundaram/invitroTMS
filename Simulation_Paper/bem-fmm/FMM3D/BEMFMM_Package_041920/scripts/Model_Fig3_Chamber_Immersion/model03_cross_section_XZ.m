%   This script computes tissue mesh intersections with the XZ plane
%   (coronal plane) located at Y=0 or elsewhere, displays and saves the
%   result
%
%   Copyright SNM 2012-2018
%   The Athinoula A. Martinos Center for Biomedical Imaging, Massachusetts General
%   Hospital & ECE Dept., Worcester Polytechnic Inst.

function model03_cross_section_XZ()
    figure();
    [~,name,~]=fileparts(pwd); L = length(name)+1;
    if ~isunix
        s = pwd; addpath(strcat(s(1:end-L), '\Engine'));
    else
        s = pwd; addpath(strcat(s(1:end-L), '/Engine'));
    end

    load CombinedMesh;

    %  Select the XZ-plane window and the y-coordinate of the XZ-plane (all in m)
    xmin = -16e-3;
    xmax = +16e-3;
    zmin = -1e-3;
    zmax = 2e-3;
    Y    = 1e-6;  % here is the plane position 
    
    %   Create coordinates of intersection contours
    tissues = length(name);
    PofXZ = cell(tissues, 1);   %   intersection nodes for a tissue
    EofXZ = cell(tissues, 1);   %   edges formed by intersection nodes for a tissue
    TofXZ = cell(tissues, 1);   %   intersected triangles
    count = [];                 %   number of every tissue present in the slice
    for m = 1:tissues
        load(name{m});   
        [P, t] = fixmesh(P, t);    
        edges = meshconnee(t);
        [Pi, ti, polymask, flag] = meshplaneintXZ(P, t, edges, Y);
        if flag % intersection found                
            count               = [count m];
            inslice             = length(count); %  total number of tissues present in the slice
            PofXZ{inslice}      = Pi;            %  intersection nodes
            EofXZ{inslice}      = polymask;      %  edges formed by intersection nodes   
            TofXZ{inslice}      = ti;            %  edges formed by intersection nodes
        end
    end

    % Display and save the result
    % Assign contour colors
    color = cool(inslice); 
    hold on;
    for m = 1:inslice
        for n = 1:size(EofXZ{m}, 1)
            i1 = EofXZ{m}(n, 1);
            i2 = EofXZ{m}(n, 2);
            line(PofXZ{m}([i1 i2], 1), PofXZ{m}([i1 i2], 3), 'Color', color(m, :), 'LineWidth', 2);
        end   
    end

    %   Insert colorbar with tissue names
    mytickmap = cell(inslice, 1);
    for m = 1:inslice
         mytickmap{m} = tissue{m};
    end
    ticks      = linspace(0.5/inslice, 1-0.5/inslice, inslice);
    colormap(color);
    colorbar;
    colorbar('Ticks', ticks, 'TickLabels', mytickmap);

    title(strcat('Coronal cross-section at y =', num2str(Y), ' m'));
    xlabel('x, m'); ylabel('z, m');
    axis 'equal';  axis 'tight';     
    axis([xmin xmax zmin zmax]);
    set(gcf,'Color','White');
    
    %  Save the result
    save('XZ_cross_section', 'PofXZ', 'EofXZ', 'TofXZ', 'count', 'color', ...
        'inslice', 'ticks', 'mytickmap', 'xmin', 'xmax', 'zmin', 'zmax', 'Y');
end

function [Pi, ti, polymask, flag] = meshplaneintXZ(P, t, edges, Y) 
%   This function implements the mesh intersection algorithm for the XZ plane 
%   Output:
%   flag = 0 -> no intersection
%   Pi - intersection nodes
%   polymask - edges formed by intersection nodes

    while 1
        index1  = P(edges(:, 1), 2)==Y; 
        index2  = P(edges(:, 2), 2)==Y; 
        indexU  = find(index1|index2);
        if ~isempty(indexU)
            Y = Y*(1+1e-12);
        else
            break;
        end 
    end    
    
    p = [0 Y 0];
    n = [0 1 0];
    flag = 1;
    Pi = [];
    ti = [];
    polymask = [];
    
    %   Find all edges (edgesI) intersecting the given plane
    index1  = P(edges(:, 1), 2)>Y; 
    index2  = P(edges(:, 2), 2)>Y;
    index3  = P(edges(:, 1), 2)<Y; 
    index4  = P(edges(:, 2), 2)<Y;
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
    %  Establish pairs of interconnected edges (pairs of interconnected edge nodes)
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

