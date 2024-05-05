function [ ] = bemf1_graphics_coil_CAD_nolight(P, t, flag) 
%   Coil 2D/3D plot with several options
%
%   Copyright SNM 2017-2019

    p = patch('vertices', P, 'faces', t);
    if flag == 0    %   non-transparent coil
        p.FaceColor = [0.4 0.94 0.94]; % [0.72 0.45 0.2];  
        %p.FaceColor = [0.5 0.5 0.5]; % [0.72 0.45 0.2];  
        p.EdgeColor = 'none';
        p.FaceAlpha = 1;
        daspect([1 1 1]);

    else
        p.FaceColor = [0.72 0.45 0.2];  
        p.EdgeColor = 'none';
        p.FaceAlpha = 0.20;      
    end
    xlabel('x, m'); ylabel('y, m'); zlabel('z, m');
    set(gcf,'Color','White');    
end