clear;clc;
%% --------------------------------generate Vi data------------------------
%------Simulation of data with real and complex components------
imdl=mk_common_model('g2s',16);% 56x56x2=6272
bkgnd= 1;
img = mk_image(imdl.fwd_model, bkgnd);

% measurement model 
load('120_i');
img.fwd_model.stimulation = stim_meas_list(meas_array);
img.fwd_model = rmfield(img.fwd_model, 'meas_select');
img.jacobian_bkgnd.value = 1;
% %% remove coarse2fine
% % save J for the first time
% q1 = img; q1fwd = rmfield(img.fwd_model,{'coarse2fine','mdl_slice_mapper'});
% q1.fwd_model = q1fwd;
% % J= calc_jacobian(calc_jacobian_bkgnd(img));
% J= calc_jacobian(calc_jacobian_bkgnd(q1));
% save J_sen J;
% 
img.solve = 'eidors_default';
vh_vol_space_1 = fwd_solve(img);
vh_vol_space = vh_vol_space_1.meas;
% figure(1);show_fem(img)
i=0;
vi_all_move = []; % save vi data
sigma_all_move = []; % save sigma data
%% 
for x = -0.75:0.0025:-0.5 %0.125  % 0.0714285
    a = x - 0.125;
    b = x + 0.125;
    c = - 0.125;
    d = 0.125;
    select_fcn = @(x,y,z)((x>a)&(x<b)&(y>c)&(y<d));  % Rectangular ((x>c)&(x<b)&(y>c)&(y<d))|((x>-0.5)&(x<-0.375)&(y>0.5)&(y<0.625))|(((x>0.5)&(x<0.625)&(y>0.5)&(y<0.625)))
    img.elem_data = 1 - 0.8 * elem_select(img.fwd_model, select_fcn);
    i=i+1;
    
    vii = fwd_solve(img);
    vi = vii.meas;
    sigma_move = img.elem_data;
    sigma_all_move = [sigma_all_move,sigma_move];
    vi_all_move = [vi_all_move,vi];
end

%% NOSER algorithm

vh_120 = vh_vol_space;
j=0;
tra_sigma = [];
for i=50:61 % sample some columns 
    vi_120 = vi_all_move(:,i);
    imb = mk_common_model('g2s', 16);
    bkgnd= 1;
    inv1d = mk_image(imb.fwd_model, bkgnd);
%   inv1d.fwd_model = imb.fwd_model;
    load('120_i');
    inv1d.jacobian_bkgnd.value= 1;
    inv1d.fwd_model.stimulation = stim_meas_list(meas_array);
    inv1d.fwd_model = rmfield(inv1d.fwd_model, 'meas_select');
    % Guass-Newton solvers
    inv1d.solve= @inv_solve_diff_GN_one_step;
    % NOSER prior
    inv1d.hyperparameter.value = .005;
    inv1d.RtR_prior=   @prior_noser;
    
    ele_double = inv1d;
    ele_double.fwd_model = rmfield(ele_double.fwd_model,{'coarse2fine','mdl_slice_mapper'});
    imgr1= inv_solve( inv1d, vh_120, vi_120);
    mm_element = imgr1.elem_data;
    tra_sigma = [tra_sigma,mm_element];    
    
    imgr1(1).calc_colours.npoints= 64;
    figure(4);
    ;
    j=j+1;subplot(3,4,j);final_center(imgr1,0.25,1);%show_slices(imgr1);
    pause(1)
end

%% save delta v

load J_sen
[a1,a2] = size(vi_all_move);
[b1,b2] = size(sigma_all_move);
delta_v = zeros(a1,a2);
for a = 1:1:a2
    delta_v(:,a) = vi_all_move(:,a) - vh_vol_space;
end
for b = 1:1:b2
    delta_sigma(:,b) = sigma_all_move(:,b) - 1;
end

save delta_vhvi_sigma-80_1fang_J delta_v J;
disp ("finished");
