
%   This computes tissue mesh intersections with the XY plane
%   (transverse plane) located at Z=0 or elsewhere, displays and saves the
%   resulting window parameters
%
%   Copyright SNM 2012-2018
%   The Athinoula A. Martinos Center for Biomedical Imaging, Massachusetts General
%   Hospital & ECE Dept., Worcester Polytechnic Inst.

function model03_cross_section_XY()

    figure();
    [~,name,~]=fileparts(pwd); L = length(name)+1;
    if ~isunix
        s = pwd; addpath(strcat(s(1:end-L), '\Engine'));
    else
        s = pwd; addpath(strcat(s(1:end-L), '/Engine'));
    end

    load CombinedMesh;
    
    %  Select the XY-plane window and the z-coordinate of the XY-plane (all in m)
    xmin = -7e-2;
    xmax = +7e-2;
    ymin = -2e-2;
    ymax = +2e-2;
    Z    =  2.0e-4; 
    
    %   Create coordinates of intersection contours and intersection edges
    tissues = length(name);
    PofXY = cell(tissues, 1);   %   intersection nodes for a tissue
    EofXY = cell(tissues, 1);   %   edges formed by intersection nodes for a tissue
    TofXY = cell(tissues, 1);   %   intersected triangles
    count = [];                 % number of every tissue present in the slice
    for m = 1:tissues
        load(name{m});   
        [P, t] = fixmesh(P, t);     
        edges = meshconnee(t);
        [Pi, ti, polymask, flag] = meshplaneintXY(P, t, edges, Z);
        if flag % intersection found                
            count               = [count m];
            inslice             = length(count); %  total number of tissues present in the slice
            PofXY{inslice}      = Pi;            %   intersection nodes
            EofXY{inslice}      = polymask;      %   edges formed by intersection nodes
            TofXY{inslice}      = ti;            %   edges formed by intersection nodes
        end
    end

    % Display contours
    color = cool(inslice);
    hold on;
    for m = 1:inslice
        for n = 1:size(EofXY{m}, 1)
            i1 = EofXY{m}(n, 1);
            i2 = EofXY{m}(n, 2);
            line(PofXY{m}([i1 i2], 1), PofXY{m}([i1 i2], 2), 'Color', color(m, :), 'LineWidth', 2);
        end   
    end

    % Insert colorbar with tissue names
    mytickmap = cell(inslice, 1);
    for m = 1:inslice
         mytickmap{m} = name{count(m)}(12:end-4);
    end
    ticks      = linspace(0.5/inslice, 1-0.5/inslice, inslice);
    colormap(color); 
    colorbar; 
    colorbar('Ticks', ticks, 'TickLabels', mytickmap);

    title(strcat('Transverse cross-section at z =', num2str(Z), ' m'));
    xlabel('x, mm'); ylabel('y, mm');
    axis 'equal';  axis 'tight';     
    axis([xmin xmax ymin ymax]);
    set(gcf,'Color','White');  
   
    box on;
    
    %   Save the result
    save('XY_cross_section', 'PofXY', 'EofXY', 'TofXY', 'count', 'color', ...
        'inslice', 'ticks', 'mytickmap', 'xmin', 'xmax', 'ymin', 'ymax', 'Z');

end

function [Pi, ti, polymask, flag] = meshplaneintXY(P, t, edges, Z) 
%   This function implements the mesh intersection algorithm for the XY plane 
%   Output:
%   flag = 0 -> no intersection
%   Pi - intersection nodes
%   polymask - edges formed by intersection nodes
%   ti - indexes into intersected triangles

    while 1
        index1  = P(edges(:, 1), 2)==Z; 
        index2  = P(edges(:, 2), 2)==Z; 
        indexU  = find(index1|index2);
        if ~isempty(indexU)
            Z = Z*(1+1e-12);
        else
            break;
        end 
    end    
    
    p = [0 0 Z];
    n = [0 0 1];
    flag = 1;
    Pi = [];
    ti = [];
    polymask = [];
    
    %   Find all edges (edgesI) intersecting the given plane
    index1  = P(edges(:, 1), 3)>Z; 
    index2  = P(edges(:, 2), 3)>Z;
    index3  = P(edges(:, 1), 3)<Z; 
    index4  = P(edges(:, 2), 3)<Z;
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

