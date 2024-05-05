function bemf1_fullmodel_graphic(COpar,P,t,Indicator)

%%% head model graphics
t0 = t(Indicator==1, :);    % skin here (change indicator if necessary)
str.EdgeColor = 'none'; str.FaceColor = [1 0.75 0.65]; str.FaceAlpha = 0.6; 
bemf2_graphics_base(P, t0, str);
t1 = t(Indicator==2, :);    % skull here (change indicator if necessary)
str1.EdgeColor = 'none'; str1.FaceColor = [1 0.75 0.65]; str1.FaceAlpha = 0.8; 
bemf2_graphics_base(P, t1, str1);

%%% Coil graphics    
Coil_All = COpar.Coil_All;

bemf1_graphics_coil_CAD(cell2mat(Coil_All(1)).P,cell2mat(Coil_All(1)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(2)).P,cell2mat(Coil_All(2)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(3)).P,cell2mat(Coil_All(3)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(4)).P,cell2mat(Coil_All(4)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(5)).P,cell2mat(Coil_All(5)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(6)).P,cell2mat(Coil_All(6)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(7)).P,cell2mat(Coil_All(7)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(8)).P,cell2mat(Coil_All(8)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(9)).P,cell2mat(Coil_All(9)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(10)).P,cell2mat(Coil_All(10)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(11)).P,cell2mat(Coil_All(11)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(12)).P,cell2mat(Coil_All(12)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(13)).P,cell2mat(Coil_All(13)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(14)).P,cell2mat(Coil_All(14)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(15)).P,cell2mat(Coil_All(15)).t, 0); hold on;
bemf1_graphics_coil_CAD_nolight(cell2mat(Coil_All(16)).P,cell2mat(Coil_All(16)).t, 0); hold on;

axis off; axis 'equal';  axis 'tight';   
daspect([1 1 1]);
set(gcf,'Color','White');
view(-180, 90);

end