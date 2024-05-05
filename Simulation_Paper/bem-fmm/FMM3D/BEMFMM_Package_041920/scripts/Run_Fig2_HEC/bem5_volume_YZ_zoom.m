%   This script accurately computes and displays electric fields sampled on
%   a cross-section (sagittal plane) via the FMM method with accurate neighbor integration
%
%   Copyright SNM/WAW 2017-2020
bem4_define_planes_zoom;
bem5_figure_scales_zoom;

%%  Load/prepare data
planeABCD = [1 0 0 X];    % Equation of the plane of the cross-section (for neighbor triangle search speedup)
%%  Post processing parameters
component   = 2;        %   field component to be plotted (1, 2, 3 or x, y, z, or 4 - total) 
temp        = ['x' 'y' 'z' 't'];
label       = temp(component);
figpos = [800 800 1600 1600];

%   sagittal plane
y = linspace(yminzoom, ymaxzoom, Ms);
z = linspace(zmin, zmax, Ms);
[Y0, Z0]  = meshgrid(y, z);
clear pointsYZ;
pointsYZ(:, 1) = X*ones(1, Ms^2);
pointsYZ(:, 2) = reshape(Y0, 1, Ms^2);
pointsYZ(:, 3) = reshape(Z0, 1, Ms^2);

%   Set up enclosing tissues (optional)
pol = 1; % 1 - skin; 2 - skull; 3 - CSF
EofYZ_closed = close_meshpolygon(EofYZ{pol}, PofYZ{pol});
poly    = meshpolygon(PofYZ{pol}, EofYZ_closed);  % cross-section in the form of an oriented polygon
in      = inpolygon(pointsYZ(:, 2), pointsYZ(:, 3), poly(:, 2), poly(:, 3));

%  Assign tissue types to observation points (required for current density plot)
obsPointTissues = assign_tissue_type_volume(pointsYZ, normals, Center, Indicator);

%% Find the E-feld at each observation point in the cross-section         
tic
Epri           = zeros(Ms*Ms, 3);
Esec           = zeros(Ms*Ms, 3);
Epri(in, :)    = bemf3_inc_field_electric(strcoil, pointsYZ(in, :), dIdt, mu0);      
Esec(in, :)    = bemf5_volume_field_electric(pointsYZ(in, :), c, P, t, Center, Area, normals, R, planeABCD);
Etotal         = Epri + Esec;   
fieldPlaneTime = toc  

%% Calculate current density at each observation point
% Observation points in free space were originally assigned tissue code 0, 
% so provide an extra "Free Space" conductivity entry and point free space obs pts to that entry.
condTemp = [cond 0];                                        
obsPointTissues(obsPointTissues == 0) = length(condTemp);

% Calculate current density J = sigma * E
condTempExpanded = transpose(condTemp(obsPointTissues));
Jtotal = repmat(condTempExpanded, 1, 3).*Etotal;

%.........................................................................
%  Plot the total E-field in the cross-section
%.........................................................................

Efig = figure('Color','w');
%  E-field plot: contour plot
if component == 4
    temp      = abs(sqrt(dot(Etotal, Etotal, 2)));
else
    temp      = (Etotal(:, component));
end
th1 = Emax;          %   in V/m
th2 = Emin;          %   in V/m
levels      = n_ELevels;
bemf2_graphics_vol_field(temp, th1, th2, levels, y, z);
xlabel('Distance y, m');
ylabel('Distance z, m');
title(['E-field (V/m), ', label, '-component in the sagittal plane.']);
 
% E-field plot: tissue boundaries
color   = prism(2); color(1, :) = [0 0 0]; color(2, :) = [0 0 0];
for m = countYZ
    edges           = EofYZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofYZ{m}(:, 2);       %   this is for the contour  
    points(:, 2)    = +PofYZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
% E-field plot: general settings 
axis 'equal';  axis 'tight';     
colormap jet; colorbar off;
axis([yminzoom ymaxzoom zmin zmax]); clim([Emin,Emax]);
grid on; set(gcf,'Color','White');
set(gcf,"Position",figpos);


%.........................................................................
%  Plot the primary E-field in the cross-section
%.........................................................................

Efig = figure('Color','w');
%  E-field plot: contour plot
if component == 4
    temp      = abs(sqrt(dot(Epri, Epri, 2)));
else
    temp      = (Epri(:, component));
end
th1 = Emax;          %   in V/m
th2 = Emin;          %   in V/m
levels      = n_ELevels;
bemf2_graphics_vol_field(temp, th1, th2, levels, y, z);
xlabel('Distance y, m');
ylabel('Distance z, m');
title(['Epri-field (V/m), ', label, '-component in the sagittal plane.']);
 
% E-field plot: tissue boundaries
color   = prism(2); color(1, :) = [0 0 0]; color(2, :) = [0 0 0];
for m = countYZ
    edges           = EofYZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofYZ{m}(:, 2);       %   this is for the contour  
    points(:, 2)    = +PofYZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
% E-field plot: general settings 
axis 'equal';  axis 'tight';     
colormap jet; colorbar off;
axis([yminzoom ymaxzoom zmin zmax]); clim([Emin,Emax]);
grid on;
set(gcf,"Position",figpos);


%.........................................................................
%  Plot the sec E-field in the cross-section
%.........................................................................

Efig = figure('Color','w');
%  E-field plot: contour plot
if component == 4
    temp      = abs(sqrt(dot(Esec, Esec, 2)));
else
    temp      = (Esec(:, component));
end
th1 = Emax;          %   in V/m
th2 = Emin;          %   in V/m
levels      = n_ELevels;
bemf2_graphics_vol_field(temp, th1, th2, levels, y, z);
xlabel('Distance y, m');
ylabel('Distance z, m');
title(['Esec-field (V/m), ', label, '-component in the sagittal plane.']);
 
% E-field plot: tissue boundaries
color   = prism(2); color(1, :) = [0 0 0]; color(2, :) = [0 0 0];
for m = countYZ
    edges           = EofYZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofYZ{m}(:, 2);       %   this is for the contour  
    points(:, 2)    = +PofYZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
% E-field plot: general settings 
axis 'equal';  axis 'tight';     
colormap jet; colorbar off;
axis([yminzoom ymaxzoom zmin zmax]); clim([Emin,Emax]);
grid on;
set(gcf,"Position",figpos);


%.........................................................................
%  Plot the vector J in the YZ cross-section
%.........................................................................

[yy,zz] = meshgrid(y,z); %Jabs = abs(sqrt(dot(Jtotal, Jtotal, 2)));

Ey = Etotal(:,2); Ez = Etotal(:,3);
Sfig = figure('Color','w'); 
imagesc(y,z,reshape(Etotal(:, component),Ms,Ms)); clim([Emin,Emax]); colormap jet;
hold on; 
L = streamslice(yy,zz,reshape(Ey,Ms,Ms),reshape(Ez,Ms,Ms),2);
set(L, 'LineWidth', 1, 'Color', 'k')
%Q = quiver(xx,yy,reshape(Jx,500,500),reshape(Jy,500,500),'Color','k','AutoScaleFactor',30);
% Display the tissue boundaries
color = ['k','k','k']; %ones(length(tissue), 3);
for m = countXY
    edges           = EofYZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofYZ{m}(:, 2);       %   this is for the contour  
    points(:, 2)    = +PofYZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m), 'LineWidth', 2.0);    %   this is contour plot
end
box on;
axis xy; axis([ymin ymax zmin zmax]);
set(gcf,"Position",figpos);
set(gca, 'Position', [0.1300 0.1100 0.7750 0.8150], 'PlotBoxAspectRatio', [3.3333 1 2.2222e+03]);
xlabel('Distance y, m');
ylabel('Distance z, m');

