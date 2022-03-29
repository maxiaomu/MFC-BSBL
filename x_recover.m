clc;clear;
load recover_x.mat
ele2 = x; 
L = size(x,2);
%%
%--------------------------------generate Vi data------------------------
%------Simulation of data with real and complex components------

img=mk_common_model('g2s',16);% 56x56x2=6272
bkgnd= 1;
imdl = mk_image(img.fwd_model, bkgnd);

% set measurement mode
load('120_i');
imdl.fwd_model.stimulation = stim_meas_list(meas_array);
imdl.fwd_model = rmfield(imdl.fwd_model, 'meas_select');
imdl.jacobian_bkgnd.value = 1;
% figure(1)

%% plot 
for i=1:L
    h1 = ele2(:,i);
    elem_ch = elem_change(h1);
    imdl.elem_data = elem_ch; 
%     final_center(imdl,0.2,0);
end
% figure(9);show_fem(imdl,[0 0 0]);
% figure(10);final_center(imdl,0.4,0);
% figure(8);show_slices(imdl);
%
figure(1);show_slices(imdl);eidors_colourbar(imdl);axis on;
% saveas(figure(1), 'f_1', 'png');
%% 三维显示
figure(2)
z=calc_slices(imdl);
h = surf(z);
shading interp;colorbar;
view(0,0);% eidors_colourbar(imdl);
% axis([0 56 0 56 0 1.2]);saveas(figure(1), 'f_2', 'png');
% axis off;
figure(3);surf(z);shading interp;colorbar;
% axis([0 56 0 56 0 1.2]);saveas(figure(3), 'f_3', 'png');
