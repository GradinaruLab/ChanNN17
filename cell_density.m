% This script detects individual cells, and calculate the minimum distance
% between neighbors and density
% Supplementary figure 6

close all 

%% img loading
%% 2nd_data
umPerPix = 0.28;

filename = 'data';
img = imread(strcat(filename, '.tif'));
figure, imshow(img(:,:,2)), title('g');

%% opening
img_open = imopen(img(:,:,2), strel('disk', 10));
figure, imshow(img_open), title('opened');

cir_range = [25, 50];
[center, radii] = imfindcircles(img_open, cir_range, 'sensitivity', 0.93);


%% remove overlapped circles
tol = min(cir_range)/2;
for i = 1:numel(radii)
    s = i+1;
    for j=s:numel(radii)
        d_ij = sqrt((center(i,1)-center(j,1)).^2+(center(i,2)-center(j,2)).^2);
        k = radii(i)+radii(j)-tol;
        
        if d_ij < k && radii(j) > 0
            if radii(i) > radii(j)
                center(j,1) = 0;
                center(j,2) = 0;
                radii(j) = 0;
            else
                center(i,1) = 0;
                center(i,2) = 0;
                radii(i) = 0;
            end
        end
    end
end
ind_cir = find(radii > 0);
center = center(ind_cir,:);
radii = radii(ind_cir);

figure, 
subplot(1,2,1), imshow(img(:,:,2)), title('original');
subplot(1,2,2), imshow(img(:,:,2)), hold on, viscircles(center, radii, 'linewidth', 1), title('opened'), hold off

dist = pdist(center);
d2 = squareform(dist);
mind = zeros(numel(radii), 1);
for i = 1:numel(mind)
    temp = d2(:,i);
    temp(find(temp==0)) = [];
    mind(i) = min(temp);
end


mind = mind.*umPerPix;

area_um = (size(img_open,1)*umPerPix) * (size(img_open,2)*umPerPix);
density = numel(mind)/area_um;
