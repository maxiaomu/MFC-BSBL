function centroids = final_center(img, amp, label)

img2 = img;
imgs = calc_slices(img2);
imgs(isnan(imgs)) = 0;

% 提取振幅大于1/4的图像
imgs_dif = max(imgs(:))-min(imgs(:));
%% 考虑在0作为中间分割点，分开导电和不导电物体
% imgs(((max(imgs(:))- amp * imgs_dif)>imgs)& (imgs>(min(imgs(:)) + amp * imgs_dif)))=0;% 去除图像中间点的像素
% imgs(imgs>(min(imgs(:)) + amp * imgs_dif))=0;% 去除图像大于某个阈值的像素
imgs(imgs<(max(imgs(:))- amp * imgs_dif))=0;% 去除图像小于某个阈值的像素
%%%%imgs(imgs>30)=0;
if isreal(imgs)
    imgs=imgs;
else
    imgs=abs(imgs);
end
% 提取图像并计算重心
[L,num] = bwlabel(imgs);
s = regionprops(L, 'centroid');
centroids = cat(1, s.Centroid);

% %%%%%%%%%%%%%%%%%%%   显示标记点和成像图
show_slices(img);
hold on;  
axis off;
%%% 有多少个中心点，个数为num
% for n=1:num
%     % 绘制重心符号：红色的“*”
%     x = centroids(n,1);
%     y = centroids(n,2);
%     plot(x,y,'r*');
%     %printf('%d',x);
%     % 在重心处标出坐标
% %     x1 = roundn(((x-65/2)/65*2),-2);  % 进行坐标变换
% %     y1 = roundn((-(y-65/2)/65*2),-2);
%     if label == 0
%         x1 = roundn(((x)/32),-2);  % 进行坐标变换
%         y1 = roundn((-(y-32)/32),-2);
%     else
% %         x1 = roundn(((x-65/2)/65*2),-2);  % 进行坐标变换
% %         y1 = roundn((-(y-65/2)/65*2),-2);
%         x1 = roundn(((x)/56),-2);  % 进行坐标变换
%         y1 = roundn((-(y-56)/56),-2);
%     end
%     str=['[',num2str(x1),',',num2str(y1),']'];
%     text(x,y,str,'fontsize',10,'Color','black');
%     hold on;
% end

%%%%%%%%%%%%%% 单独现实标记点
% for n=1:num
%     % 绘制重心符号：红色的“*”
%     x = centroids(n,1);
%     y = centroids(n,2);
%     
%     %printf('%d',x);
%     % 在重心处标出坐标
%     x1 = roundn((x-65/2)/65*2,-2)+1;  % 进行坐标变换
%     y1 = roundn(-(y-65/2)/65*2,-2);
%     plot(x1,y1,'r.','MarkerSize',16);
%     axis([0,1,0,1])
%     str=['[',num2str(x1),',',num2str(y1),']'];
%     text(x1,y1,str,'fontsize',10,'Color','black');
%     hold on;
% end