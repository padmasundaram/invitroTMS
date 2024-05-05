function [targ_dirs, as] = target_def(cursor_info,n_targ_dirs,P_tmp,t_tmp)

targ_dirs=zeros(1,1,3); 
    
t_tmp_an=[1 0 0]; 
n_tmp=cursor_info.Position;
ex_coil_in_mri=((t_tmp_an)-dot(t_tmp_an,-n_tmp)*-n_tmp)/(norm((t_tmp_an)-dot(t_tmp_an,-n_tmp)*-n_tmp));

theta=(pi/(180/n_targ_dirs));
%R=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
R_rot=rotationmat3D(theta,n_tmp);
   
%%%ROTATING WITH A SPECIFIED AMOUNT
targ_dirs(1,1,:)=(R_rot*ex_coil_in_mri')';
%targ_dirs = targ_dirs/norm(targ_dirs);   
 
figure()
patch('Faces',t_tmp,'Vertices',P_tmp,'FaceColor','y','Edgecolor','none','FaceAlpha',0.2);hold on; 
axis off;camlight; lighting gouraud; daspect([ 1 1 1]); 
set(gcf,'Color','White'); colormap jet; view(-180,90);
quiver3(cursor_info.Position(1),cursor_info.Position(2),cursor_info.Position(3),squeeze(targ_dirs(1,:,1))',squeeze(targ_dirs(1,:,2))',squeeze(targ_dirs(1,:,3))',0.1);

[as bs] = knnsearch(P_tmp,cursor_info.Position,'K',100,'Distance','euclidean');%% creating the target vetcor for MNE optimization
%quiver3(P_tmp(as,1),P_tmp(as,2),P_tmp(as,3),squeeze(targ_dirs(1,1))',squeeze(targ_dirs(1,2))',squeeze(targ_dirs(1,3))',0.1);
end