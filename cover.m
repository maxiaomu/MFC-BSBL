clc;clear;

load recover_x.mat
% x(x<0.01)=0;
% x((x>0.01)&(x<0.1))=0.8;
% x(x>0.01)=1;
ele2 = x; % 适用于真实的数据
% ele2 = x/2; % 适用于仿真无噪声的数据
L = size(x,2);
%%
%--------------------------------生成正向Vi数据------------------------
%------Simulation of data with real and complex components------
% img=mk_common_model('e2s',16);% 32*32
img=mk_common_model('g2s',16);% 56x56x2=6272
bkgnd= 1;
imdl = mk_image(img.fwd_model, bkgnd);

%设置激励测量模式
load('120_i');
imdl.fwd_model.stimulation = stim_meas_list(meas_array);
imdl.fwd_model = rmfield(imdl.fwd_model, 'meas_select');
imdl.jacobian_bkgnd.value = 1;
% figure(1)
% for i=1:9
%     img.elem_data = ele1(:,i); 
%     show_slices(img);picname=['ori_' num2str(i) '.fig'];saveas(gcf,picname);%,[inf,inf,0,1,1]
% end
% figure(2)
% for i=1:9
%     img.elem_data = ele2(:,i); 
%     show_slices(img);picname=['rec_' num2str(i) '.fig'];saveas(gcf,picname);%,[inf,inf,0,1,1] 
% end
% %% 绘制在一张图
% for i=1:9
%     img.elem_data = ele1(:,i); 
%     figure(1);subplot(3,3,i);show_slices(img);
% end
%% 绘制正向模型 
for i=1:L
    h1 = ele2(:,i);
    elem_ch = elem_change(h1);
    imdl.elem_data = elem_ch; 
%     img.elem_data = ele2(:,i); 
%     figure(4);subplot(2,4,i); 
    final_center(imdl,0.2,0);
%     show_slices(imdl);
%     final_center(imdl,0.2,0);
%     show_fem(imdl,[0 0 0]);
end
% figure(9);show_fem(imdl,[0 0 0]);
% figure(10);final_center(imdl,0.4,0);
% figure(8);show_slices(imdl);
%
figure(10);show_slices(imdl);eidors_colourbar(imdl);axis on;
% saveas(figure(10), 'dandian2_51', 'png');
%% 三维显示
figure(13)
z=calc_slices(imdl);
% c=calc_colours(z,imdl);
h = surf(z);
% set(h,'edgecolor','none');
shading interp;colorbar;
view(0,0);% eidors_colourbar(imdl);
% axis([0 56 0 56 0 1.2]);saveas(figure(13), 'dandian2_52', 'png');
% axis off;
figure(14);surf(z);shading interp;colorbar;
% axis([0 56 0 56 0 1.2]);saveas(figure(14), 'dandian2_53', 'png');
%% 保存ele_data数据
% ele_50 = abs(imdl.elem_data);
% save ele_50 ele_50
%% NOSER prior 
% imdl.jacobian_bkgnd.value= 1;
% imdl.fwd_model.stimulation = stim_meas_list( meas_array );
% imdl.fwd_model = rmfield( imdl.fwd_model, 'meas_select');  
%  
% load deltav_single6;
% delta_v = deltav';  % delta_v 120*10
% vh=delta_v(:,10)-delta_v(:,10);
% imdl.solve= @inv_solve_diff_GN_one_step;
% imdl.hyperparameter.value = .03;
% imdl.RtR_prior=@prior_noser;%prior_laplace,prior_gaussian_HPF,
% q1 = imdl;
% q1.fwd_model = rmfield(q1.fwd_model,{'coarse2fine','mdl_slice_mapper'});
% for i=1:10
%     figure(4);subplot(3,4,i)
%     vi=delta_v(:,i);
%     imgr=inv_solve(q1, vh, vi);
%     imgr.elem_data = ele2(:,i);
%     show_slices(imgr);%final_center(imgr,0.25,1); 
% end
% % load deltav_single7
% % deltav = deltav';
% % for i=1:1:10
% %     figure(5);subplot(3,4,i)
% %     vi=deltav(:,i);
% %     imgr=inv_solve(q1, vh, vi);
% %     show_slices(imgr);%final_center(imgr,0.25,1); 
% % end