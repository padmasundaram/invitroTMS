%   This script accurately computes and displays electric fields sampled on
%   a cross-section (transverse plane) via the FMM method with accurate neighbor integration
%
%   Copyright SNM/WAW 2017-2020

%%  Load/prepare data
planeABCD = [0 0 1 Z];        % Equation of the plane of the cross-section (Ax + By + Cz + D = 0)(meters) for neighbor triangle search acceleration

%%  Post processing parameters
component   = 2;        %   field component to be plotted (1, 2, 3 or x, y, z, or 4 - total) 
temp        = ['x' 'y' 'z' 't'];
label       = temp(component);


%   Coronal plane
x = linspace(xmin, xmax, Ms);
y = linspace(ymin, ymax, Ms);
[X0, Y0]  = meshgrid(x, y);
% clear pointsXY;
pointsXY       = zeros(Ms*Ms,3);
pointsXY(:, 1) = reshape(X0, 1, Ms^2);
pointsXY(:, 2) = reshape(Y0, 1, Ms^2);  
pointsXY(:, 3) = Z*ones(1, Ms^2);

%  Set up enclosing tissues (optional: suppresses visualization of saturated E-field outside the head model)
pol = 1; % 1 - skin; 2 - skull; 3 - CSF
EofXY_closed = close_meshpolygon(EofXY{pol}, PofXY{pol});
poly    = meshpolygon(PofXY{pol}, EofXY_closed);  % cross-section in the form of an oriented polygon
in      = inpolygon(pointsXY(:, 1), pointsXY(:, 2), poly(:, 1), poly(:, 2));

%  Assign tissue types to observation points (required for current density plot)
obsPointTissues = assign_tissue_type_volume(pointsXY, normals, Center, Indicator);

%%  Find the E-field at each observation point in the cross-section        
tic
%pointsXY       = 1e-3*pointsXY;     % Convert back to m
Epri           = zeros(Ms*Ms, 3);
Esec           = zeros(Ms*Ms, 3);
Epri(in, :)    = bemf3_inc_field_electric(strcoil, pointsXY(in, :), dIdt, mu0);      
Esec(in, :)    = bemf5_volume_field_electric(pointsXY(in, :), c, P, t, Center, Area, normals, R, planeABCD);
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

Efig = figure;

% E-field plot: contour plot
if component == 4
    temp      = abs(sqrt(dot(Etotal, Etotal, 2)));
else
    temp      = (Etotal(:, component));
end
th1 = Emax;          %   in V/m
th2 = Emin;          %   in V/m
levels      = n_ELevels;
bemf2_graphics_vol_field(temp, th1, th2, levels, x, y);
xlabel('Distance x, m');
ylabel('Distance y, m');
title(['E-field (V/m), ', label, '-component in the XY plane.']);

% E-field plot: tissue boundaries
color   = prism(length(tissue)); color(2, :) = [1 1 1];
for m = countXY
    edges           = EofXY{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXY{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXY{m}(:, 2);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
% E-field plot: general settings 
axis 'equal';  axis 'tight';     
colormap jet; colorbar; clim([Emin,Emax]);
axis([xmin xmax ymin ymax]);
grid on; set(gcf,'Color','White');
set(gcf,"Position",[81 431 1014 1017]);
EtotXY = temp;

%.........................................................................
%  Plot the primary E-field in the cross-section
%.........................................................................

Efig = figure;

% E-field plot: contour plot
if component == 4
    temp      = abs(sqrt(dot(Etotal, Etotal, 2)));
else
    temp      = (Epri(:, component));
end
th1 = Eprimax;          %   in V/m
th2 = Eprimin;          %   in V/m
levels      = n_ELevels;
bemf2_graphics_vol_field(temp, th1, th2, levels, x, y);
xlabel('Distance x, m');
ylabel('Distance y, m');
title(['Epri-field (V/m), ', label, '-component in the XY plane.']);

% E-field plot: tissue boundaries
color   = prism(length(tissue)); color(2, :) = [1 1 1];
for m = countXY
    edges           = EofXY{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXY{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXY{m}(:, 2);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
% E-field plot: general settings 
axis 'equal';  axis 'tight';     
colormap jet; colorbar; clim([Eprimin,Eprimax]);
axis([xmin xmax ymin ymax]);
grid on; set(gcf,'Color','White');
set(gcf,"Position",[81 431 1014 1017]);
EpriXY = temp;

%.........................................................................
%  Plot the sec E-field in the cross-section
%.........................................................................

Efig = figure;
% E-field plot: contour plot
if component == 4
    temp      = abs(sqrt(dot(Etotal, Etotal, 2)));
else
    temp      = (Esec(:, component));
end
th1 = Emax;          %   in V/m
th2 = Emin;          %   in V/m
levels      = n_ELevels;
bemf2_graphics_vol_field(temp, th1, th2, levels, x, y);
xlabel('Distance x, m');
ylabel('Distance y, m');
title(['Esec-field (V/m), ', label, '-component in the XY plane.']);

% E-field plot: tissue boundaries
color   = prism(length(tissue)); color(2, :) = [1 1 1];
for m = countXY
    edges           = EofXY{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXY{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXY{m}(:, 2);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
% E-field plot: general settings 
axis 'equal';  axis 'tight';     
colormap jet; colorbar; clim([Emin,Emax]);
axis([xmin xmax ymin ymax]);
grid on; set(gcf,'Color','White');
set(gcf,"Position",[81 431 1014 1017]);
EsecXY = temp;


%.........................................................................
%  Plot the total, pri and sec E-field profiles in the XY cross-section
%.........................................................................

EtotXYsq = reshape(EtotXY, Ms, Ms); 
EpriXYsq = reshape(EpriXY, Ms, Ms); 
EsecXYsq = reshape(EsecXY, Ms, Ms); 

NX = Ms/2; 
figure('color','w'); 
L = plot(y, EtotXYsq(:,NX), 'r', y, EpriXYsq(:,NX), 'g', y, EsecXYsq(:,NX), 'b');

set(L(1), 'LineWidth', 2, 'Color', [228 26 28]/255); % red
set(L(2), 'LineWidth', 2, 'Color', [77 175 74]/255); % green
set(L(3), 'LineWidth', 2, 'Color', [55 126 184]/255);% blue

hold on;
L1 = line([-0.002 -0.002], [-45 45]); set(L1, 'LineStyle', '--', 'Color', 'k');
L2 = line([+0.002 +0.002], [-45 45]); set(L2, 'LineStyle', '--', 'Color', 'k');
xlabel('y (m)');
ylabel('E (V/m)');

ylim([-45,45]);

