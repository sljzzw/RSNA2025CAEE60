function out = CSsim(Img, CSfactor, varargin)
% CSsim - Simulate Compressed Sensing MRI Reconstruction (Cartesian)
%
%   out = CSsim(Img, CSfactor, 'Name', value, ...)
%
%   This function performs compressed sensing reconstruction using
%   variable density random sampling and TV + L1 regularization.
%
%   Inputs:
%       Img        - Input image (2D complex image)
%       CSfactor   - Undersampling factor (e.g. 5 => 20% k-space)
%
%   Optional Name-Value Pairs:
%       'MaxOuterIter'    - Outer CG iterations (default = 4)
%       'PolyPower'       - Polynomial power for PDF (default = 20)
%       'TVWeight'        - TV regularization weight (default = 0.1)
%       'XFMWeight'       - L1 regularization weight (default = 0.1)
%       'Itnlim'          - Inner CG iterations (default = 30)
%       'Show'            - Show figures (default = 1)
%       'SaveIntermediate' - Save intermediate steps (default = true)
%
%   Outputs:
%       out.img_full, out.img_lr, out.img_dc, out.img_cs
%       out.mask, out.pdf, out.k_full, out.k_us
%       out.MSE, out.L1, out.actpctg
%
%   Author: Zhongwei Zhang, MD, PhD
%   Washington University School of Medicine
%   2025

%% --- Parse Inputs ---
p = inputParser;
addOptional(p, 'MaxOuterIter', 100, @isnumeric);
addOptional(p, 'PolyPower', 20, @isnumeric);
addOptional(p, 'TVWeight', 0.1, @isnumeric);
addOptional(p, 'XFMWeight', 0.1, @isnumeric);
addOptional(p, 'Itnlim', 30, @isnumeric);
addOptional(p, 'Show', 1, @isnumeric);
addOptional(p, 'SaveIntermediate', true, @islogical);
parse(p, varargin{:});
opt = p.Results;

DN = size(Img);
pctg = 1 / CSfactor;

%% --- Sampling PDF and Mask ---
pdf = genPDF(DN, opt.PolyPower, pctg, 2, 0, 0);
[mask, ~, actpctg] = genSampling(pdf, 500, 2);
FT = p2DFT(mask, DN, 1, 2);

%% --- K-space Generation ---
% Full-grid k-space
k_full_grid = fft2c(Img);    

% Sampled data (vector)
data = FT * Img;

% Reconstruct undersampled k-space grid from sampled vector
if ~isequal(size(mask), DN)
    error('mask size %s must match image size %s.', mat2str(size(mask)), mat2str(DN));
end

% p2DFT returned full-length vector, reshape directly
k_us_grid = reshape(data, DN) .* mask;




k_full = k_full_grid;
k_us   = k_us_grid;

%% --- Transform Operator ---
XFM = 1;
TVWeight = opt.TVWeight;
xfmWeight = opt.XFMWeight;

%% --- Initialization ---
im_dc = FT' * (data ./ pdf);
im_init = im_dc;

%% --- Visualization Setup ---
if opt.Show
    set(0,'DefaultFigureColor','k'); % black background
    figure(100); clf;
    subplot(2,3,1), imshow(abs(Img),[]), title('Original','Color','w');
    subplot(2,3,2), imshow(log(abs(fftshift(k_full))+1),[]), title('Full k-space (log)','Color','w');
    subplot(2,3,3), imshow(log(abs(fftshift(k_us))+1),[]), title('Undersampled k-space (log)','Color','w');
    subplot(2,3,4), imshow(mask,[]), title('Sampling Mask','Color','w');
    subplot(2,3,5), imshow(abs(im_dc),[]), title('Zero-fill DC','Color','w');
    drawnow;
end

%% --- Reconstruction Parameters ---
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight = TVWeight;
param.xfmWeight = xfmWeight;
param.Itnlim = opt.Itnlim;

res = XFM * im_init;
MSE = [];
L1 = [];
intermediate = {};

%% --- Iterative Reconstruction ---
for n = 1:opt.MaxOuterIter
    res = fnlCg(res, param);
    param.TVWeight = param.TVWeight * 0.8;
    param.xfmWeight = param.xfmWeight * 0.8;

    im_res = XFM' * res;
    tmpmse = sum(abs(param.FT * im_res - data).^2,'all');
    tmpl1 = sum(abs(param.TV * im_res),'all') + sum(abs(res),'all');

    MSE(end+1) = tmpmse;
    L1(end+1) = tmpl1;

    if opt.SaveIntermediate
        intermediate{n} = im_res;
    end

    if opt.Show
        subplot(2,3,6), imshow(abs(im_res),[]); 
        title(sprintf('CS Iter %d', n),'Color','w');
        figure(101); clf;
        plot(MSE,'r','LineWidth',1.5); hold on;
        plot(L1,'c','LineWidth',1.5);
        legend('MSE','L1+TV','TextColor','w');
        set(gca,'Color','k','XColor','w','YColor','w');
        title('Cost Function','Color','w');
        drawnow;
    end
end

%% --- Low-resolution reference ---
mask_lr = genLRSampling_pctg(DN, pctg, 1, 0);
im_lr = ifft2c(zpad(fft2c(Img).*mask_lr, DN(1), DN(2)));
im_full = ifft2c(zpad(fft2c(Img), DN(1), DN(2)));

if opt.Show
    figure(200); clf;
    imshow(abs(cat(2,im_full,im_lr,im_dc,im_res)),[]);
    title(sprintf('Full | Low Res | Zero-fill DC | CS (%.1f%% sampling)', pctg*100),'Color','w');
end


%% --- Outputs ---
out.img_full = im_full;
out.img_lr = im_lr;
out.img_dc = im_dc;
out.img_cs = im_res;
out.mask = mask;
out.pdf = pdf;
out.k_full = k_full;
out.k_us = k_us;
out.MSE = MSE;
out.L1 = L1;
out.actpctg = actpctg;

if opt.SaveIntermediate
    out.intermediate = intermediate;
end

% ==========================
% Visualization from "out"
% ==========================
img_full = abs(out.img_full);
img_dc   = abs(out.img_dc);
img_cs   = abs(out.img_cs);
mask     = out.mask;
k_full   = out.k_full;
k_us     = out.k_us;

% --- Normalize for display ---
img_full = img_full ./ max(img_full(:));
img_dc   = img_dc ./ max(img_dc(:));
img_cs   = img_cs ./ max(img_cs(:));

k_full_mag = log(1 + abs(k_full));
k_full_mag = k_full_mag ./ max(k_full_mag(:));

k_us_mag = log(1 + abs(k_us));
k_us_mag = k_us_mag ./ max(k_us_mag(:));

mask_disp = double(mask);

% --- Summary Figure ---
figure('Color','k');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

% Row 1: K-space
nexttile;
imshow(k_full_mag, []);
title('Full K-space','Color','w','FontSize',14);
axis off; colormap gray;

nexttile;
imshow(k_us_mag, []);
title('Sparse K-space','Color','w','FontSize',14);
axis off; colormap gray;

nexttile;
imshow(mask_disp, []);
title('Sampling Mask','Color','w','FontSize',14);
axis off; colormap gray;

% Row 2: Images
nexttile;
imshow(img_full, []);
title('Full Image','Color','w','FontSize',14);
axis off; colormap gray;

nexttile;
imshow(img_dc, []);
title('Direct Recon (ZF+DC)','Color','w','FontSize',14);
axis off; colormap gray;

nexttile;
imshow(img_cs, []);
title('CS Recon Image','Color','w','FontSize',14);
axis off; colormap gray;

sgtitle(sprintf('CS Recon Simulation (CS factor = %g, pctg = %.1f%%)', ...
    CSfactor, out.actpctg*100), 'Color','y','FontSize',16);

end


%{

Plot W-space image

% Input: img_cs (CS reconstructed image)
img_cs = abs(img_cs);
img_cs = img_cs ./ max(img_cs(:));

wname   = 'db4';
nLevels = 3;

% --- Wavelet decomposition ---
[C, S] = wavedec2(img_cs, nLevels, wname);
[LL, LH, HL, HH] = dwt2(img_cs, wname);

% --- Original W-space visualization ---
Wspace = [LL, LH; HL, HH];
Wspace = abs(Wspace);
Wspace = log(1 + Wspace);
Wspace = Wspace ./ max(Wspace(:));

% --- Denoising in W-space ---
sigma = median(abs(C)) / 0.6745;
T = sigma * sqrt(2*log(numel(C)));
C_denoised = wthresh(C, 's', T);
img_denoised = waverec2(C_denoised, S, wname);

% --- Denoised W-space visualization ---
[LLd, LHd, HLd, HHd] = dwt2(img_denoised, wname);
Wspace_denoised = [LLd, LHd; HLd, HHd];
Wspace_denoised = abs(Wspace_denoised);
Wspace_denoised = log(1 + Wspace_denoised);
Wspace_denoised = Wspace_denoised ./ max(Wspace_denoised(:));

% --- K-space of denoised image ---
K_denoised = fftshift(fft2(img_denoised));
K_denoised_mag = log(1 + abs(K_denoised));
K_denoised_mag = K_denoised_mag ./ max(K_denoised_mag(:));

% === Final Figure ===
figure('Color','k');
tiledlayout(1,5,'Padding','compact','TileSpacing','compact');

nexttile;
imshow(img_cs, []);
title('Original','Color','w','FontSize',14);
axis off; colormap gray;

nexttile;
imshow(Wspace, []);
title('W-space','Color','w','FontSize',14);
axis off; colormap gray;

nexttile;
imshow(Wspace_denoised, []);
title('Denoised W-space','Color','w','FontSize',14);
axis off; colormap gray;

nexttile;
imshow(img_denoised, []);
title('Denoised Image','Color','w','FontSize',14);
axis off; colormap gray;

nexttile;
imshow(K_denoised_mag, []);
title('K-space (Denoised)','Color','w','FontSize',14);
axis off; colormap gray;

sgtitle('Wavelet Domain and K-space Visualization','Color','w','FontSize',16);



%}





function mask = genLRSampling_pctg(imSize,pctg,distType,disp)
% mask = genLRSampling_pctg(imSize,pctg,distType,disp)
%
% Creates a circular mask around the center
% user enters percentage of pixels out of the total number.
%
% (c) Michael Lustig 2007


sx = imSize(1);
sy = imSize(2);


if sum(imSize==1)==0  % 2D	
	[x,y] = meshgrid(linspace(-1,1,sy),linspace(-1,1,sx));
	switch distType
		case 1
			r = max(abs(x),abs(y));
		otherwise
			r = sqrt(x.^2+y.^2);
			r = r/max(abs(r(:)));			
	end
else %1d
	r = abs(linspace(-1,1,max(sx,sy)));
end

[nothing, circOrder] = sort(r(:));

mask = zeros(imSize);
mask(circOrder(1:floor(pctg*sx*sy))) = 1;

if disp
	figure, 
	subplot(211), imshow(mask)
	if sum(imSize==1)==0
		subplot(212), plot(mask(end/2+1,:));
	else
		subplot(212), plot(mask);
	end
end

end



function res = ifft2c(x)
%
%
% res = ifft2c(x)
% 
% orthonormal centered 2D ifft
%
% (c) Michael Lustig 2005

res = sqrt(length(x(:)))*ifftshift(ifft2(fftshift(x)));


end


function res = ifft2c(x)
%
%
% res = ifft2c(x)
% 
% orthonormal centered 2D ifft
%
% (c) Michael Lustig 2005

res = sqrt(length(x(:)))*ifftshift(ifft2(fftshift(x)));

end

