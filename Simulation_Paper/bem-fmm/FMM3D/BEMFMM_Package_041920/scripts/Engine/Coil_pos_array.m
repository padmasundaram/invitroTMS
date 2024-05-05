function [Coil_arr] = Coil_pos_array(lx,ly)

xu = [1 0 0];yu = [0 1 0];zu = [0 0 1];

r_s = 0.095;a = ly + 2; b = lx - 2;
XY = [-a a;-a b;-a -b;-a -a;-b a; -b b;-b -b;-b -a;b a;b b;b -b;b -a;a a;a b;a -b;a -a]*1e-2;

Z = sqrt(r_s^2-sum((XY.^2),2));
MovXYZ(:,1:2) = XY;
MovXYZ(:,3) = Z;
for i = 1:16
ang= [atan2(norm(cross(MovXYZ(i,:),xu)),dot(MovXYZ(i,:),xu)) atan2(norm(cross(MovXYZ(i,:),yu)),dot(MovXYZ(i,:),yu)) atan2(norm(cross(MovXYZ(i,:),zu)),dot(MovXYZ(i,:),zu))]
Nxyz(i,:) = cos(ang);
end
Coil_arr.MovXYZ = MovXYZ;
Coil_arr.Nxyz = Nxyz;
end