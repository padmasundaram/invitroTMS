
% choice of tissue_dim: 4x4, 4x2, 4x1 (for Fig5)
% For supplementary:
% 1x1, 1x2, 1x4
% 2x1, 2x2, 2x4
% 4x1, 4x2, 4x4
% just set tissue_dim accordingly
%

tissue_dim = '4x4'; %<---- 

bem0_load_model;
bem1_setup_coil;
bem2_charge_engine;
bem4_define_planes;
bem5_figure_scales;
bem5_volume_XY;
bem5_volume_YZ;

