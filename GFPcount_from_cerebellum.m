% This script detects calbindin-positive cells and NLS-GFP from cerebellum images to calculate %GFP and mean intensity
% Figure 2, Supplementary figure 2

close all 
clear all

%% img loading
img = imread('data.tif');

% figure, 
% subplot(1,3,1), imshow(img(:,:,1)), title('r');
% subplot(1,3,2), imshow(img(:,:,2)), title('g');
% subplot(1,3,3), imshow(img), title('merged');
% saveas(gcf, strcat(foldername, '1_original.jpg'));


%% parameters
gfp_thresh = .05;
areath = .1;


%% opening
img_open = imopen(img(:,:,1), strel('disk', 5));
figure,
subplot(1,2,1), imshow(img(:,:,1)), title('r');
subplot(1,2,2), imshow(img_open), title('r-opened');

[center, radii] = imfindcircles(img_open, [15, 35], 'sensitivity', 0.87);
figure, imshow(img(:,:,1)), hold on, viscircles(center, radii, 'linewidth', 1), title('opened'), hold off

mask_red = zeros(size(img(:,:,1)));
calbindin_ispos = zeros(numel(radii), 1);
g_mean_inten = zeros(numel(radii), 1);
g_bg_inten = zeros(numel(radii), 1);
green = img(:,:,2);
for i=1:numel(radii)
    temp = zeros(size(img(:,:,1)));
    temp(round(center(i,2)), round(center(i,1))) = 1;
    cir = strel('disk', floor(radii(i)));
    conv_cir = conv2(temp, double(cir.Neighborhood), 'same');
    mask_red = mask_red | conv_cir;
    ind_cir = find(conv_cir > 0);
        
    % gfp (nucleus) detection    
    mask_sq = ones(size(cir.Neighborhood));
    conv_sq = conv2(temp, mask_sq, 'same');
    ind_sq = find(conv_sq > 0);
    
    if numel(ind_sq) ~= numel(mask_sq)
        g_mean_inten(i) = nan;
        calbindin_ispos(i) = nan;
        continue
    end
    
    g_calpos = reshape(conv_cir(ind_sq), size(mask_sq));
    g_fltrd_sq = reshape(green(ind_sq), size(mask_sq));
    
    % background subtraction for each cell area in GFP image
    g_bg = imgaussfilt(g_fltrd_sq, size(g_fltrd_sq, 1)*5);
    g_bg = g_bg - (min(g_bg(:))-min(g_fltrd_sq(:)));
    g_bg_inten(i) = double(max(g_bg(:)))./254;
    g_fltrd_sq2 = g_fltrd_sq - g_bg;
    
    % nucleus detection after bg subtraction
    g_fltrd_bw = imbinarize(g_fltrd_sq2 , 'global');
    g_fltrd = g_fltrd_bw & cir.Neighborhood;
    g_fltrd_ch = bwconvhull(g_fltrd, 'object', 4);
    
    se_dia = strel('diamond',1);
    g_fltrd_sm = imerode(g_fltrd_ch, se_dia);
    g_fltrd_sm = imerode(g_fltrd_sm, se_dia);
    
    cc = bwconncomp(g_fltrd_sm, 4); 
    g_fltrd_max = zeros(size(g_fltrd_sm));
    
    if cc.NumObjects >= 1
        idx = 1;
        if cc.NumObjects > 1
            num_pix = cellfun(@numel,cc.PixelIdxList);
            [biggest,idx] = max(num_pix);
        end
        g_fltrd_max(cc.PixelIdxList{idx}) = 1; 
    end
        
    g_fltrd_dl = imdilate(g_fltrd_max, se_dia);
    g_fltrd_dl = imdilate(g_fltrd_dl, se_dia);
    
    g_fltrd_areath = bwareaopen(g_fltrd_dl, round(pi*radii(i)^2*areath), 4);

    g_outline = bwperim(g_fltrd_areath);
    g_fltrd_bdry = g_fltrd_sq2;
    g_fltrd_bdry(g_outline) = 255;

    ind_inside = find(g_fltrd_areath > 0);
    g_mean_inten(i) = median(g_fltrd_sq2(ind_inside));
    
    if isnan(g_mean_inten(i))
        g_mean_inten(i) = 0;
    end
    
%         figure,
%         subplot(1,3,1), imagesc(double(g_fltrd_sq)./254), title('green'), ...
%             set(gca, 'yticklabel', []), set(gca, 'xticklabel', []), ...
%             colorbar, axis square ij;
%         subplot(1,3,2), imagesc(double(g_bg)./254), title('gauss - bg'), ...
%             set(gca, 'yticklabel', []), set(gca, 'xticklabel', []), ...
%             colorbar, axis square ij;
%         subplot(1,3,3), imagesc(double(g_fltrd_sq2)./254), title('bg subtracted'), ...
%             set(gca, 'yticklabel', []), set(gca, 'xticklabel', []), ...
%             colorbar, axis square ij;
%         saveas(gcf, strcat(foldername, sprintf('4_bgsubtraction_%d.jpg', i)));
%     
%         figure, 
%         subplot(2,3,1), imshow(g_fltrd_sq2), title(sprintf('%d', i))
%         subplot(2,3,2), imshow(g_fltrd), title('ostu');
%         subplot(2,3,3), imshow(g_fltrd_ch), title('convex hull');
%         subplot(2,3,4), imshow(g_fltrd_sm), title('erode')
%         subplot(2,3,5), imshow(g_fltrd_max), title('areath');
%         subplot(2,3,6), imshow(g_fltrd_bdry), title(sprintf('mean: %.2f', g_mean_inten(i)/254));
%         saveas(gcf, strcat(foldername, sprintf('5_gfp_%d.jpg', i)));
%         pause
 end


calbindin_ispos(isnan(calbindin_ispos)) = [];
g_mean_inten(isnan(g_mean_inten)) = [];

g_mean_inten = g_mean_inten./254;

gfp_mean_inten = mean(g_mean_inten(find(g_mean_inten >= gfp_thresh)));


disp('')
disp('---- Result')
fprintf(' # of calbindin+ = %d\n', numel(calbindin_ispos));
fprintf(' # of GFP+       = %d\n', numel(find(g_mean_inten >= gfp_thresh)));
fprintf(' GFP/calbindin   = %f\n', numel(find(g_mean_inten >= gfp_thresh))/numel(calbindin_ispos)*100);
fprintf(' GFP intensity   = %f\n', gfp_mean_inten);