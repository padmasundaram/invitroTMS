function [Coil_repos] = position_coils(Coil_init)

strcoil_All = [];
Coil_All = [];

coilaxis        = [0 0 1];
theta = 0;
for ij =1:Coil_init.len
    
Coil_init.strcoil_tmp.Pwire   = meshrotate2(Coil_init.strcoil.Pwire, coilaxis, theta);
Coil_init.Coil_tmp.P          = meshrotate2(Coil_init.Coil.P, coilaxis, theta);
%  Tilt the coil axis itself as required
Nx = Coil_init.Nxyz(ij,1); 
Ny = Coil_init.Nxyz(ij,2); 
Nz = Coil_init.Nxyz(ij,3);

Coil_init.strcoil_tmp.Pwire = meshrotate1(Coil_init.strcoil_tmp.Pwire, Nx, Ny, Nz);
Coil_init.Coil_tmp.P        = meshrotate1(Coil_init.Coil_tmp.P, Nx, Ny, Nz);
%  Move the coil as required
MoveX = Coil_init.MovXYZ(ij,1); 
MoveY = Coil_init.MovXYZ(ij,2); 
MoveZ = Coil_init.MovXYZ(ij,3);

Coil_init.strcoil_tmp.Pwire = Coil_init.strcoil_tmp.Pwire + [MoveX MoveY MoveZ];

Coil_init.Coil_tmp.P   = Coil_init.Coil_tmp.P + [MoveX MoveY MoveZ];
%   Define the observation line from the bottom center of the coil into the
%   head (for field plotting)
M = 10000;        
argline      = linspace(0, 100e-3, M);              %   distance along a 100 mm long line   
dirline      = -[Nx Ny Nz]/norm([Nx Ny Nz]);        %   line direction (along the coil axis)   
offline      = 0e-3;                                %   offset from the coil
pointsline_tmp(1:M, 1) = MoveX + dirline(1)*(argline + offline);
pointsline_tmp(1:M, 2) = MoveY + dirline(2)*(argline + offline);
pointsline_tmp(1:M, 3) = MoveZ + dirline(3)*(argline + offline);

%%% saving the results
Coil_repos.strcoil_All{ij} = Coil_init.strcoil_tmp;
Coil_repos.Coil_All{ij} = Coil_init.Coil_tmp;
end
end