% This script detects individual cells and calculate the number of colors
% expressed in each cell.
% Figure 4f

close all 
clear all

%% img loading
img = imread('data.tif');
img_raw = img;
r_open = 7;
cir_range = [10,20];
sens = 0.86; 
foldername = 'data/';
thresh = 0.2;


figure, 
subplot(2,2,1), imshow(img(:,:,1)), title('r');
subplot(2,2,2), imshow(img(:,:,2)), title('g');
subplot(2,2,3), imshow(img(:,:,3)), title('b');
subplot(2,2,4), imshow(img), title('merged');


%% contrast adjustment for each channel
img_ad = zeros(size(img), 'uint16');
for i=1:3
    img_ad(:,:,i) = imadjust(img(:,:,i));
end


figure,
subplot(2,2,1), imshow(img_ad(:,:,1)), title('r-ad');
subplot(2,2,2), imshow(img_ad(:,:,2)), title('g-ad');
subplot(2,2,3), imshow(img_ad(:,:,3)), title('b-ad');
subplot(2,2,4), imshow(img_ad), title('merged-ad');


%% rgb 2 hsv
img_hsv = rgb2hsv(img_ad);
figure,
subplot(2,2,1), imshow(img_hsv(:,:,1)), title('h');
subplot(2,2,2), imshow(img_hsv(:,:,2)), title('s');
subplot(2,2,3), imshow(img_hsv(:,:,3)), title('v');
subplot(2,2,4), imshow(img_hsv), title('merged');

img_gray = img_hsv(:,:,3);

%% opening
img_open = imopen(img_gray, strel('disk', r_open));
figure,
subplot(1,2,1), imshow(img_gray), title('gray');
subplot(1,2,2), imshow(img_open), title('gray-opened');


%% CHT for each channel
center = [];
radii =[];
for i=1:3
    [c_ch, r_ch] = imfindcircles(imopen(img_ad(:,:,i), strel('disk', r_open)), cir_range, 'sensitivity', sens);
    figure,
    subplot(1,2,1), imshow(img_ad(:,:,i))
    subplot(1,2,2), imshow(img_ad(:,:,i)), hold on, viscircles(c_ch, r_ch, 'linewidth', 1), hold off;
    center = [center; c_ch];
    radii = [radii; r_ch];
end

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
subplot(1,2,1), imshow(img_ad), title('img_ad');
subplot(1,2,2), imshow(img_ad), title('cell detected'), hold on, 
for i=1:numel(radii)
    viscircles(center(i,:), radii(i)-1, 'linewidth', 1),
    text(center(i,1), center(i,2), sprintf('%d', i), 'color', 'r');
end
hold off


%% normalize intensity to the range of 0-1
img_norm = zeros(size(img));
for i = 1:3
    temp = double(img_ad(:,:,i));
    norm = temp - min(temp(:));
    norm = norm./max(norm(:));
    img_norm(:,:,i) = norm;
end

figure,
subplot(1,4,1), imshow(img_norm(:,:,1)), title('r-norm'), colorbar;
subplot(1,4,2), imshow(img_norm(:,:,2)), title('g-norm'), colorbar;
subplot(1,4,3), imshow(img_norm(:,:,3)), title('b-norm'), colorbar;
subplot(1,4,4), imshow(img_norm), title('merged');
figure, imshow(img_norm), title('norm')


img_open = zeros(size(img));
for i = 1:3
    img_open(:,:,i) = imopen(img_norm(:,:,i), strel('disk', r_open));
end

figure,
subplot(2,2,1), imshow(img_open(:,:,1)), title('r-norm'), colorbar;
subplot(2,2,2), imshow(img_open(:,:,2)), title('g-norm'), colorbar;
subplot(2,2,3), imshow(img_open(:,:,3)), title('b-norm'), colorbar;
subplot(2,2,4), imshow(img_open), title('merged');

figure, imshow(img_open), title('opened_cell-detected'), 


hoverlapped = zeros(numel(radii), 1);
meandist = zeros(numel(radii), 3);
distance = zeros(numel(radii), 3);
binedge = 0:0.01:1;
cell_inten_total = [];
for i=1:numel(radii)
    temp = zeros(size(img(:,:,1)));
    temp(round(center(i,2)), round(center(i,1))) = 1;
    cir = strel('disk', floor(radii(i)-1));
    conv_cir = conv2(temp, double(cir.Neighborhood), 'same');
    
    ind_cir = find(conv_cir > 0);
    
    cell_inten = zeros(numel(ind_cir), 3);
    for j=1:3
        norm_1ch = img_open(:,:,j);
        cell_inten(:, j) = norm_1ch(ind_cir);
    end
     
    pdfit = zeros(numel(binedge), 3);
    for j=1:3
        pd = fitdist(cell_inten(:,j), 'normal');
        pdfit(:,j) = pdf(pd, binedge);
        meandist(i,j) = pd.mu;
        pd.mu
        mean(cell_inten(:,j))
    end
    
    if i <= 10
        figure,
        histogram(cell_inten(:,1), binedge, 'normalization', 'probability', 'displaystyle', 'stairs', 'edgecolor', 'r'), ...
            title(sprintf('cell: %d', i)); hold on;
        histogram(cell_inten(:,2), binedge, 'normalization', 'probability', 'displaystyle', 'stairs', 'edgecolor', 'g');
        histogram(cell_inten(:,3), binedge, 'normalization', 'probability', 'displaystyle', 'stairs', 'edgecolor', 'b');
        plot(binedge, pdfit(:,1)./sum(pdfit(:,1)), 'r', 'linewidth', 2);
        plot(binedge, pdfit(:,2)./sum(pdfit(:,2)), 'g', 'linewidth', 2);
        plot(binedge, pdfit(:,3)./sum(pdfit(:,3)), 'b', 'linewidth', 2);
        line([thresh, thresh], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k');
        hold off
     end
    
    
     xlswrite(strcat(foldername, sprintf('cell_inten_%d.xlsx', i)), cell_inten, 1, sprintf('B%d', 3*(i-1)+1));

end


%% counting the number of expressed channels
exnum = zeros(numel(radii), 1);
for i=1:numel(radii)
    exnum(i) = numel(find(meandist(i,:)>=thresh));
end

result = [numel(find(exnum==1)), numel(find(exnum==2)), numel(find(exnum==3))];