function [ ] = bemf1_graphics_coil_CAD_3color(coil, flag) 
    Gc = [0.11 0.93 0.1];
    Rc = [1 0 0];
    Bc = [0.01 0.01 0.98];
    if flag == 1    %   x coil
        Pfull = coil.P;
        p = patch('vertices', coil.PX, 'faces', coil.tX);
        p.FaceColor = Rc; % [0.72 0.45 0.2];  
        p.EdgeColor = 'none';
        p.FaceAlpha = 1;
        daspect([1 1 1]);
        camlight; 
        lighting flat;   
    elseif flag == 2
        p = patch('vertices', coil.PY, 'faces', coil.tY);
        p.FaceColor = Gc; % [0.72 0.45 0.2];  
        p.EdgeColor = 'none';
        p.FaceAlpha = 1;
        daspect([1 1 1]);

    elseif flag == 3
        p = patch('vertices', coil.PZ, 'faces', coil.tZ);
        p.FaceColor = Bc; % [0.72 0.45 0.2];  
        p.EdgeColor = 'none';
        p.FaceAlpha = 1;
        daspect([1 1 1]);
%         camlight; 
%         lighting flat; 
    end
    xlabel('x, m'); ylabel('y, m'); zlabel('z, m');
    set(gcf,'Color','White');    
end