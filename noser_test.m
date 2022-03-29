clc;clear;
%% 
load simulation_data/delta_vhvi_sigma-80_1fang_J.mat

% Y = Phi_ * X_ /4.0;
vi_120 = delta_v(:,80);
vh_0 = vi_120-vi_120;
%% 
% load vhvi_sigma90_1fang.mat
% vi_120 = vi_all_move(:,51);
% vh_0 = vh_vol_space;
%%
imdl=mk_common_model('g2s',16);% 56x56x2=6272
bkgnd= 1;
img = mk_image(imdl.fwd_model, bkgnd);

load('120_i');
img.fwd_model.stimulation = stim_meas_list(meas_array);
img.fwd_model = rmfield(img.fwd_model, 'meas_select');
img.jacobian_bkgnd.value = 1;

vh_00 = fwd_solve(img);
vh = vh_00;
vh.meas = vh_0;
vi = vh_00;
vi.meas = vi_120;

% Create Inverse Model
inv2d= eidors_obj('inv_model', 'EIT inverse');
inv2d.reconst_type= 'difference';
inv2d.jacobian_bkgnd.value= 1;
% This is not an inverse crime; inv_mdl != fwd_mdl
inv2d.fwd_model= img.fwd_model;
%%
inv2d.solve= @inv_solve_diff_GN_one_step;
inv2d.hyperparameter.value = .003;
inv2d.RtR_prior= @prior_noser;
%%    
ele_double = img;
ele_double.fwd_model = rmfield(ele_double.fwd_model,{'coarse2fine','mdl_slice_mapper'});
imgr1= inv_solve(inv2d, vh, vi);
figure(1);show_slices(imgr1);eidors_colourbar(imgr1);
