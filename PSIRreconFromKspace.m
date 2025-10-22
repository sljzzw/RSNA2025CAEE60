function psir_img = PSIRreconFromKspace(IR_img, Ref_img, varargin)
% PSIR_RECON_MULTIFRAME - Phase-Sensitive Inversion Recovery Reconstruction
%   for multi-frame or cine IR MRI datasets
%
%   psir_img = psir_recon_multiframe(IR_img, Ref_img)
%   psir_img = psir_recon_multiframe(IR_img, Ref_img, 'normalize', true)
%
% DESCRIPTION:
%   This function performs phase-sensitive inversion recovery (PSIR)
%   reconstruction across multiple frames. For each frame, the sign of the
%   magnetization is recovered using the reference phase. It preserves
%   tissue polarity and improves contrast compared to magnitude IR.
%
% INPUTS:
%   IR_img   : complex IR dataset
%              Size can be [Nx Ny Nslice Nframe] or [Nx Ny Nframe]
%
%   Ref_img  : complex reference dataset
%              Can be:
%                 1) Single reference image (same size as [Nx Ny Nslice])
%                 2) Per-frame reference (same size as IR_img)
%
% OPTIONAL PARAMETERS:
%   'normalize' : logical (default = true)
%                 Normalize output to [-1, 1] per frame for display.
%
% OUTPUT:
%   psir_img : reconstructed PSIR image
%              Same size as IR_img
%
% METHOD:
%   For each frame:
%     1. Compute magnitude of IR image
%     2. Compute phase sign using reference phase
%     3. Multiply magnitude by sign
%     4. Optional per-frame normalization
%
% REFERENCE:
%   Kellman P, Arai AE. “Phase-sensitive inversion-recovery for detecting 
%   myocardial infarction using gadolinium-delayed hyperenhancement.”
%   Magn Reson Med. 2001;45(5):846–852.
%
% AUTHOR:
%   Zhongwei Zhang, MD, PhD
%   Washington University School of Medicine
%   Email: zhongweiz@wustl.edu
%{


% S_IR and S_REF are coil-combined complex images
phi_IR  = angle(S_IR);
phi_REF = angle(S_REF);

% Phase difference
dphi = phi_IR - phi_REF;
dphi = angle(exp(1i*dphi));  % wrap to [-pi, pi]

% Sign determination
sign_map = sign(cos(dphi));
sign_map(sign_map==0) = 1;

% PSIR image
PSIR = abs(S_IR) .* sign_map;

% Optional scaling
PSIR_scaled = PSIR ./ abs(S_REF) * 1000;



% Suppose IR_img and Ref_img are complex arrays [Nx Ny Nslice Nframe]
load IR_complex.mat  % example: IR_img
load Ref_complex.mat % example: Ref_img

% PSIR reconstruction across all frames
psir_out = psir_recon_multiframe(IR_img, Ref_img, 'normalize', true);

% Visualize one slice across time
z = 1;
figure;
for f = 1:size(psir_out,4)
    subplot(3,ceil(size(psir_out,4)/3),f);
    imshow(psir_out(:,:,z,f), [-1 1]);
    title(['Frame ' num2str(f)]);
end
colormap(gray);
sgtitle('PSIR Reconstruction - Multi-frame');


%}

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
addParameter(p, 'normalize', true, @islogical);
parse(p, varargin{:});
doNormalize = p.Results.normalize;

% -------------------------------------------------------------------------
% Check dimensions
% -------------------------------------------------------------------------
IR_size = size(IR_img);
ndim = ndims(IR_img);

if ndim == 2
    % Single frame 2D image
    IR_img = reshape(IR_img, IR_size(1), IR_size(2), 1, 1);
elseif ndim == 3
    % Could be Nx×Ny×Nframe or Nx×Ny×Nz
    if size(Ref_img,3) == 1 || size(Ref_img,3) == size(IR_img,3)
        IR_img = reshape(IR_img, IR_size(1), IR_size(2), 1, IR_size(3));
    else
        IR_img = reshape(IR_img, IR_size(1), IR_size(2), IR_size(3), 1);
    end
elseif ndim ~= 4
    error('IR_img must be 2D, 3D or 4D.');
end

[Nx, Ny, Nz, Nf] = size(IR_img);

% Match Ref_img dimensions
if ismatrix(Ref_img)
    Ref_img = reshape(Ref_img, Nx, Ny, 1, 1);
elseif ndims(Ref_img) == 3
    Ref_img = reshape(Ref_img, Nx, Ny, Nz, 1);
elseif ndims(Ref_img) == 4 && ~isequal(size(Ref_img), size(IR_img))
    error('Ref_img size must match IR_img or be single reference.');
end

% -------------------------------------------------------------------------
% Initialize output
% -------------------------------------------------------------------------
psir_img = zeros(Nx, Ny, Nz, Nf, 'like', IR_img);

% -------------------------------------------------------------------------
% Frame loop
% -------------------------------------------------------------------------
for f = 1:Nf
    for z = 1:Nz
        % IR and reference for this frame and slice
        IR_frame  = IR_img(:,:,z,f);
        if size(Ref_img,4) > 1
            Ref_frame = Ref_img(:,:,z,f);
        else
            Ref_frame = Ref_img(:,:,z,1);
        end

        % Magnitude
        mag_IR = abs(IR_frame);

        % Sign from phase projection
        sign_val = sign(real(IR_frame .* conj(Ref_frame)));
        sign_val(sign_val == 0) = 1;

        % PSIR reconstruction
        psir_frame = mag_IR .* sign_val;

        % Optional normalization
        if doNormalize
            max_val = max(psir_frame(:));
            if max_val ~= 0
                psir_frame = psir_frame ./ max_val;
            end
        end

        % Store
        psir_img(:,:,z,f) = psir_frame;
    end
end

% If original input was 2D or 3D, reshape back
if ndim == 2
    psir_img = squeeze(psir_img(:,:,1,1));
elseif ndim == 3 && size(IR_size,2) == 3
    if size(Ref_img,4) == 1
        psir_img = squeeze(psir_img(:,:,1,:));
    end
end

end
