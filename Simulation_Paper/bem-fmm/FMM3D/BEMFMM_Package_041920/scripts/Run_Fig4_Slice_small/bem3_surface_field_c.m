%   This script plots the induced surface charge density for
%   any brain compartment surface (plots the density + optionally
%   coil geometry).
%
%   Copyright SNM/WAW 2017-2020


w = 5e-9;

figure;
subplot(121);

%%   Graphics
tissue_to_plot = 'Ringer';

objectnumber    = find(strcmp(tissue, tissue_to_plot));
temp            = eps0*c(Indicator==objectnumber);  % the real charge density is eps0*c

bemf2_graphics_surf_field(P, t, temp, Indicator, objectnumber);
title(strcat('Solution: Surface charge density in C/m^2 for: ', tissue{objectnumber}));

% % Coil centerline graphics 
% bemf1_graphics_coil_CAD(Coil.P, Coil.t, 1);
% hold on;
% plot3(pointsline(:, 1), pointsline(:, 2), pointsline(:, 3), '-m', 'lineWidth', 5);

% General
axis tight; view(30, 35), camzoom(1);
clim([-w,w]);

brighten(0.4);

% (fnz == 0) && (fnx > 0)
% (fnz == 0) && (fnx < 0)
fnA = {};
cnt = 1;
for ii = 1:length(Indicator)
    if (Indicator(ii) == objectnumber)
        if dot(normals(ii,:),[0,1,0]) == 1
            fnA{cnt} = ii; cnt = cnt + 1;
        end
    end
end
fnA = cell2mat(fnA);

c_fnA = eps0*c(fnA);

fnB = {};
cnt = 1;
for ii = 1:length(Indicator)
    if (Indicator(ii) == objectnumber)
        if dot(normals(ii,:),[0,-1,0]) == 1
            fnB{cnt} = ii; cnt = cnt + 1;
        end
    end
end
fnB = cell2mat(fnB);

c_fnB = eps0*c(fnB);

save('chargedensity_ringer.mat', 'c_fnA', 'c_fnB', 'P', 't', 'temp', 'Indicator', 'objectnumber');


subplot(122);
%%   Graphics
tissue_to_plot = 'Slice';

objectnumber    = find(strcmp(tissue, tissue_to_plot));
temp            = eps0*c(Indicator==objectnumber);  % the real charge density is eps0*c
bemf2_graphics_surf_field(P, t, temp, Indicator, objectnumber);
title(strcat('Solution: Surface charge density in C/m^2 for: ', tissue{objectnumber}));

% % Coil centerline graphics 
% bemf1_graphics_coil_CAD(Coil.P, Coil.t, 1);
% hold on;
% plot3(pointsline(:, 1), pointsline(:, 2), pointsline(:, 3), '-m', 'lineWidth', 5);

% General
axis tight; view(30, 35), camzoom(1);

brighten(0.4);
clim([-w,w]);

% (fnz == 0) && (fnx > 0)
% (fnz == 0) && (fnx < 0)
fnA = {};
cnt = 1;
for ii = 1:length(Indicator)
    if (Indicator(ii) == objectnumber)
        if dot(normals(ii,:),[0,1,0]) == 1
            fnA{cnt} = ii; cnt = cnt + 1;
        end
    end
end
fnA = cell2mat(fnA);

c_fnA = eps0*c(fnA);

fnB = {};
cnt = 1;
for ii = 1:length(Indicator)
    if (Indicator(ii) == objectnumber)
        if dot(normals(ii,:),[0,-1,0]) == 1
            fnB{cnt} = ii; cnt = cnt + 1;
        end
    end
end
fnB = cell2mat(fnB);

c_fnB = eps0*c(fnB);

save('chargedensity_slice.mat', 'c_fnA', 'c_fnB', 'P', 't', 'temp', 'Indicator', 'objectnumber');
