function centroids = final_center(img, amp, label)

img2 = img;
imgs = calc_slices(img2);
imgs(isnan(imgs)) = 0;

% Remove irrelevant data points
imgs_dif = max(imgs(:))-min(imgs(:));
%% 
% imgs(((max(imgs(:))- amp * imgs_dif)>imgs)& (imgs>(min(imgs(:)) + amp * imgs_dif)))=0;% remove 
% imgs(imgs>(min(imgs(:)) + amp * imgs_dif))=0;
imgs(imgs<(max(imgs(:))- amp * imgs_dif))=0;
%%%%imgs(imgs>30)=0;
if isreal(imgs)
    imgs=imgs;
else
    imgs=abs(imgs);
end
% Calculate the center of gravity
[L,num] = bwlabel(imgs);
s = regionprops(L, 'centroid');
centroids = cat(1, s.Centroid);

% %%%%%%%%%%%%%%%%%%%   figure
show_slices(img);
hold on;  
axis off;
%%% number of point
% for n=1:num
%     % red“*”
%     x = centroids(n,1);
%     y = centroids(n,2);
%     plot(x,y,'r*');
%     %printf('%d',x);
%     % mark
%     if label == 0
%         x1 = roundn(((x)/32),-2);  % 进行坐标变换
%         y1 = roundn((-(y-32)/32),-2);
%     else
%         x1 = roundn(((x)/56),-2);  % 进行坐标变换
%         y1 = roundn((-(y-56)/56),-2);
%     end
%     str=['[',num2str(x1),',',num2str(y1),']'];
%     text(x,y,str,'fontsize',10,'Color','black');
%     hold on;
% end
