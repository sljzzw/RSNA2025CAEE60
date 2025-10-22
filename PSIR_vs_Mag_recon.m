% ----------------------------------------------------------------------------------
% | Tissue              | Pre-contrast T1 (ms) | Post-contrast T1 (ms, ~10–20 min) |
% | ------------------- | -------------------: | --------------------------------: |
% | Blood               |                 1500 |                           300–450 |
% | Myocardium (normal) |             950–1000 |                           350–500 |
% | Scar (fibrosis/LGE) |            1100–1400 |                           200–350 |
% | CSF                 |                 4000 |                         3000–4000 |
% ---------------------------------------------------------------------------------
% Inversion Recovery Gradient Echo Simulation
% Simulate inversion recovery signal evolution at 1.5T for different tissues
% before and after Gadolinium contrast
%
% Signal model:
%   S(t) = M0 * (1 - 2 * exp(-t ./ T1))
%
% Reference T1 values at 1.5T:
%   Myocardium (native): ~950 ms  (Messroghli DR, MRM 2003)
%   Blood (native):      ~1500 ms (Gai ND, JMRI 2017)
%   Scar (post MI):      ~450 ms  (Kellman P, JCMR 2007)
%   CSF:                 ~4000 ms (Stanisz GJ, MRM 2005)
%
% Post-contrast (typical 10–15 min after GBCA):
%   Myocardium: ~450 ms
%   Blood:      ~300 ms
%   Scar:       ~250 ms
%   CSF:        ~2500 ms  (minimal enhancement)
%
% TI = 300 ms (typical late gadolinium enhancement)
%
% Author: Zhongwei Zhang, MD, PhD
% 
%% Inversion Recovery GRE Signal Simulation (Mag vs PSIR)
% Simulate IR-GRE signal at TI = 300 ms for different tissues
% Before and after gadolinium contrast enhancement.
% Reference:
%   Kim RJ et al. Circulation. 1999;100:1992–2002.
%   Kellman P, Arai AE. J Cardiovasc Magn Reson. 2007;9(3):525–537.
%   Simonetti OP et al. Radiology. 2001;218(1):215–223.

clear; clc;

%% Parameters
t = 0:1:800;             % time in ms
TI = 300;                % inversion time in ms
M0 = 1;                  % equilibrium magnetization

% Pre-contrast T1 (ms) at 1.5T
T1_blood_pre       = 1200;
T1_myocardium_pre  = 950;
T1_scar_pre        = 500;
T1_csf_pre         = 2500;

% Post-contrast T1 (ms)
T1_blood_post       = 350;
T1_myocardium_post  = 450;
T1_scar_post        = 250;
T1_csf_post         = 1800;

%% Signal model: IR-GRE
S_pre.blood       = M0 * (1 - 2*exp(-t/T1_blood_pre));
S_pre.myocardium  = M0 * (1 - 2*exp(-t/T1_myocardium_pre));
S_pre.scar        = M0 * (1 - 2*exp(-t/T1_scar_pre));
S_pre.csf         = M0 * (1 - 2*exp(-t/T1_csf_pre));

S_post.blood       = M0 * (1 - 2*exp(-t/T1_blood_post));
S_post.myocardium  = M0 * (1 - 2*exp(-t/T1_myocardium_post));
S_post.scar        = M0 * (1 - 2*exp(-t/T1_scar_post));
S_post.csf         = M0 * (1 - 2*exp(-t/T1_csf_post));

%% Plot colors
col_blood = [0 0.5 0];          % dark green
col_myocardium = [0.8 0.6 0];   % dark yellow
col_scar = [1 0 0];             % red
col_csf = [0 0 1];              % blue

figure('Color','k');

%% Subplot 1 - Without Contrast
subplot(2,1,1)
hold on; grid on;

% --- Plot magnitude reconstruction (solid)
plot(t, abs(S_pre.blood), 'Color', col_blood, 'LineWidth', 2);
plot(t, abs(S_pre.myocardium), 'Color', col_myocardium, 'LineWidth', 2);
plot(t, abs(S_pre.scar), 'Color', col_scar, 'LineWidth', 2);
plot(t, abs(S_pre.csf), 'Color', col_csf, 'LineWidth', 2);

% --- Plot PSIR reconstruction (dashed)
plot(t, S_pre.blood, '--', 'Color', col_blood, 'LineWidth', 1.5);
plot(t, S_pre.myocardium, '--', 'Color', col_myocardium, 'LineWidth', 1.5);
plot(t, S_pre.scar, '--', 'Color', col_scar, 'LineWidth', 1.5);
plot(t, S_pre.csf, '--', 'Color', col_csf, 'LineWidth', 1.5);

% TI and zero lines
xline(TI,'--w','LineWidth',1.5);
yline(0,'w-','LineWidth',1);

set(gca,'Color','k','XColor','w','YColor','w','FontSize',14);
xlabel('Time (ms)','Color','w');
ylabel('Signal Intensity','Color','w');
title('Without Contrast','Color','w');
xlim([0 800]); ylim([-1.2 1.2]);

%% Subplot 2 - With Contrast
subplot(2,1,2)
hold on; grid on;

% --- Mag Recon (solid)
plot(t, abs(S_post.blood), 'Color', col_blood, 'LineWidth', 2);
plot(t, abs(S_post.myocardium), 'Color', col_myocardium, 'LineWidth', 2);
plot(t, abs(S_post.scar), 'Color', col_scar, 'LineWidth', 2);
plot(t, abs(S_post.csf), 'Color', col_csf, 'LineWidth', 2);

% --- PSIR Recon (dashed)
plot(t, S_post.blood, '--', 'Color', col_blood, 'LineWidth', 1.5);
plot(t, S_post.myocardium, '--', 'Color', col_myocardium, 'LineWidth', 1.5);
plot(t, S_post.scar, '--', 'Color', col_scar, 'LineWidth', 1.5);
plot(t, S_post.csf, '--', 'Color', col_csf, 'LineWidth', 1.5);

% TI and zero lines
xline(TI,'--w','LineWidth',1.5);
yline(0,'w-','LineWidth',1);

set(gca,'Color','k','XColor','w','YColor','w','FontSize',14);
xlabel('Time (ms)','Color','w');
ylabel('Signal Intensity','Color','w');
title('With Gadolinium Contrast','Color','w');
xlim([0 800]); ylim([-1.2 1.2]);

%% Custom legend with dummy handles
subplot(2,1,1); % put legend on top subplot
hold on;

% Line style legend (Mag vs PSIR)
h_mag = plot(nan, nan, 'w-', 'LineWidth', 2);
h_psir = plot(nan, nan, 'w--', 'LineWidth', 1.5);

% Color legend (tissues)
h_blood = plot(nan, nan, '-', 'Color', col_blood, 'LineWidth', 2);
h_myo = plot(nan, nan, '-', 'Color', col_myocardium, 'LineWidth', 2);
h_scar = plot(nan, nan, '-', 'Color', col_scar, 'LineWidth', 2);
h_csf = plot(nan, nan, '-', 'Color', col_csf, 'LineWidth', 2);

legend([h_mag h_psir h_blood h_myo h_scar h_csf], ...
       {'Mag Recon','PSIR Recon','Blood','Myocardium','Scar','CSF'}, ...
       'TextColor','w','Location','southeast');
