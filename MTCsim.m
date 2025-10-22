function out = MTCsim(varargin)
% MTCsim - MT simulation with offset-frequency sweep, multi-tissue, and MTR maps
% -------------------------------------------------------------------------
% Usage:
%   out = MTCsim(Name,Value,...)
%
% Important Name-Value pairs (defaults shown):
%   'OffsetHz'    : single offset frequency used for detailed Mz plots (3000)
%   'OffsetSweep' : vector of offsets [Hz] to sweep and compute MTR (e.g. 500:250:5000)
%   'PulseDur'    : MT pulse duration (s) default 0.012
%   'B1uT'        : peak B1 amplitude in microTesla (default 15)
%   'PulseShape'  : 'gauss' (only) 
%   'PulseSamples': samples across pulse (default 5000)
%   'Tissues'     : cell array of tissue names (default {'Blood','Myocardium','Scar'})
%   'UserTissue'  : struct to add as 'user' tissue
%   'dt'          : integrator step (s) default pulseDur/pulseSamples
%   'alpha'       : imaging flip deg used for MTR (default 30 deg)
%   'PlotLogX'    : true/false for log-x MTR plot (default false)
%   'RecoveryTime': time after MT pulse before readout (s) (default 0)
%   'ShowPlots'   : true/false to display plots (default true)
%   'ShowLineshape': true/false to show bound pool lineshape (default false)
%
% Output:
%   out.MTR_vs_offset (Ntissue x Noffset)
%   out.Mz_free_curve, out.Mz_bound_curve (for selected/first offset)
%   out.Mxy_free_curve (transverse magnetization magnitude)
%   out.* (other diagnostic fields)
%
% Author: Zhongwei Zhang, MD, PhD
% Date: 2025-10-12 (improved)

% --- parse inputs (vendor-informed pragmatic defaults) -------------------
p = inputParser;
% Use 1.2 kHz as a pragmatic Siemens/GE-style product default (see Cohen-Adad et al. / MPM guidance).
addParameter(p,'OffsetHz',2000,@(x)validateattributes(x,{'numeric'},{'scalar','finite','nonnan'}));
addParameter(p,'OffsetSweep',[],@(x)validateattributes(x,{'numeric'},{'vector','finite'}));
% PulseDur set to 10 ms (typical MTprep/MTsat pulses in the literature are 5–30 ms; 10–12 ms is common).
addParameter(p,'PulseDur',15e-3,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
% Peak B1 in microTesla: use a moderate default (2–5 µT are typical starting points to reach strong saturation).
addParameter(p,'B1uT',12,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
addParameter(p,'PulseShape','gauss',@(x)ischar(x) || isstring(x));
addParameter(p,'PulseSamples',5000,@(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',3}));
addParameter(p,'Tissues',{'Blood','Myocardium','Scar'});
addParameter(p,'UserTissue',struct(),@isstruct);
addParameter(p,'alpha',10,@(x)validateattributes(x,{'numeric'},{'scalar'}));
addParameter(p,'PlotLogX',false,@islogical);
addParameter(p,'RecoveryTime',0,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
addParameter(p,'ShowPlots',true,@islogical);
addParameter(p,'ShowLineshape',false,@islogical); % New: lineshape visualization
parse(p,varargin{:});
par = p.Results;


% Constants
gamma = 2*pi*42.577e6; % rad/s/T

%% --- build normalized pulse envelope (peak=1) and RF waveforms
[env, tp, dt, B1_t, omega1_t, flip_eff_deg, B1_area] = buildPulseEnvelope(par, gamma);

% Create full time vector including initial state (N+1 samples)
tp_full = [tp(1)-dt, tp];    % initial time just before first sample
Ntp = numel(tp_full);

%% --- tissue definitions
TissueDB = buildTissueDB(par.Tissues, par.UserTissue);
TissueDB = verifyTissueParameters(TissueDB); % Auto-corrects parameters
Ntissue = numel(TissueDB);

%% --- prepare offsets to simulate
[offsets, Noffsets, selected_offset, selected_idx] = setupOffsets(par);

% Preallocate outputs (time length N+1)
MTR_vs_offset = zeros(Ntissue, Noffsets);
Mz_free_curve_sel = zeros(Ntp, Ntissue);
Mz_bound_curve_sel = zeros(Ntp, Ntissue);
Mxy_free_curve_sel = zeros(Ntp, Ntissue);
Mx_free_curve_sel = zeros(Ntp, Ntissue);
My_free_curve_sel = zeros(Ntp, Ntissue);

%% --- Compute bound pool lineshapes if requested
if par.ShowLineshape
    [lineshape_absorption, lineshape_dispersion, linewidths_Hz] = computeLineshapes(TissueDB, offsets);
else
    lineshape_absorption = [];
    lineshape_dispersion = [];
    linewidths_Hz = [];
end

%% --- main offset sweep loop
for io = 1:Noffsets
    offsetHz = offsets(io);
    Delta = 2*pi*offsetHz; % rad/s (off-resonance)
    
    for ti = 1:Ntissue
        params_t = TissueDB(ti);
        
        % initial magnetization at equilibrium (normalized M0f + M0b = 1)
        M0 = [0; 0; (1 - params_t.M0b); 0; 0; params_t.M0b];
        
        % Integrate through MT pulse - returns transient curves for selected offset
        [M_final, Mz_free_temp, Mz_bound_temp, Mx_free_temp, My_free_temp, Mxy_free_temp] = ...
            integrateBlochMcConnell(M0, omega1_t, Delta, params_t, dt, tp, io==selected_idx);
        
        % If recovery time specified, integrate further with no RF
        if par.RecoveryTime > 0
            M_final = applyRecovery(M_final, Delta, params_t, dt, par.RecoveryTime);
        end
        
        % Save curves for selected offset
        if io == selected_idx
            Mz_free_curve_sel(:,ti) = Mz_free_temp;
            Mz_bound_curve_sel(:,ti) = Mz_bound_temp;
            Mx_free_curve_sel(:,ti) = Mx_free_temp;
            My_free_curve_sel(:,ti) = My_free_temp;
            Mxy_free_curve_sel(:,ti) = Mxy_free_temp;
        end

        % Compute MTR
        MTR_vs_offset(ti, io) = computeMTR(M_final, params_t, par.alpha);
    end
end

%% --- package outputs
out = struct();
out.params = par;
out.flip_eff_deg = flip_eff_deg;
out.B1_area = B1_area;
out.offsets = offsets;
out.MTR_vs_offset = MTR_vs_offset;
out.Mz_free_curve = Mz_free_curve_sel;
out.Mz_bound_curve = Mz_bound_curve_sel;
out.Mx_free_curve = Mx_free_curve_sel;
out.My_free_curve = My_free_curve_sel;
out.Mxy_free_curve = Mxy_free_curve_sel;
out.t = tp_full;           % pulse time vector (s) including initial sample
out.B1_t = [0; B1_t];      % pad B1_t with 0 at t0 to match N+1 (no RF at t0)
out.omega1_t = [0; omega1_t];
out.TissueDB = TissueDB;
out.selected_offset = selected_offset; 
out.selected_idx = selected_idx;
out.lineshape_absorption = lineshape_absorption;
out.lineshape_dispersion = lineshape_dispersion;
out.linewidths_Hz = linewidths_Hz;

%% --- plotting
if par.ShowPlots
    createSummaryPlots(out, par, selected_offset);
end

%% --- print brief summary
printSimulationSummary(out, TissueDB);

end

%% Helper Functions

function [env, tp, dt, B1_t, omega1_t, flip_eff_deg, B1_area] = buildPulseEnvelope(par, gamma)
% Build normalized pulse envelope and compute RF parameters
    N = par.PulseSamples;
    tp = linspace(-par.PulseDur/2, par.PulseDur/2, N);
    
    switch lower(char(par.PulseShape))
        case 'gauss'
            % choose sigma so that ±(PulseDur/2) ~ 3*sigma -> sigma = PulseDur/6
            sigma = par.PulseDur/6;
            env = exp(-(tp.^2)/(2*sigma.^2));
        otherwise
            error('Unsupported pulse shape: %s', par.PulseShape);
    end
    env = env(:);
    env = env / max(env); % normalized to 1
    dt = par.PulseDur / (N-1);
    
    % B1(t) in Tesla (peak = par.B1uT microTesla)
    B1_t = (par.B1uT * 1e-6) * env;               % T
    omega1_t = gamma * B1_t;                      % rad/s

    % compute B1 area and effective flip
    B1_area = trapz(tp, B1_t);                    % T·s
    flip_eff_rad = gamma * B1_area;
    flip_eff_deg = flip_eff_rad * 180/pi;
end

function TissueDB = buildTissueDB(names, user_tissue)
% Build tissue parameter database
    if ischar(names) || isstring(names)
        names = {char(names)};
    end
    TissueDB = struct([]);
    cnt = 0;
    
    for ii = 1:numel(names)
        nm = lower(char(names{ii}));
        cnt = cnt + 1;
        switch nm
            case 'blood'
                TissueDB(cnt).Name = 'Blood';
                TissueDB(cnt).M0b = 0.03;      % Bound pool fraction
                TissueDB(cnt).T1f = 1.8;       % Free pool T1 (s)
                TissueDB(cnt).T2f = 0.18;      % Free pool T2 (s)
                TissueDB(cnt).T1b = 1.0;       % Bound pool T1 (s)
                TissueDB(cnt).T2b = 1e-5;      % Bound pool T2 (s) - very short
                TissueDB(cnt).kf = 1.5;        % Forward exchange rate (s^-1)
                TissueDB(cnt).kb = TissueDB(cnt).kf * TissueDB(cnt).M0b / (1 - TissueDB(cnt).M0b);
                
            case 'myocardium'
                TissueDB(cnt).Name = 'Myocardium';
                TissueDB(cnt).M0b = 0.10;
                TissueDB(cnt).T1f = 1.3;
                TissueDB(cnt).T2f = 0.052;
                TissueDB(cnt).T1b = 1.0;
                TissueDB(cnt).T2b = 1e-5;
                TissueDB(cnt).kf = 3.0;
                TissueDB(cnt).kb = TissueDB(cnt).kf * TissueDB(cnt).M0b / (1 - TissueDB(cnt).M0b);
                
            case {'scar','fibrotic','fibrosis'}
                TissueDB(cnt).Name = 'Fibrotic Scar';
                TissueDB(cnt).M0b = 0.15;
                TissueDB(cnt).T1f = 1.2;
                TissueDB(cnt).T2f = 0.05;
                TissueDB(cnt).T1b = 1.0;
                TissueDB(cnt).T2b = 1e-5;
                TissueDB(cnt).kf = 4.0;
                TissueDB(cnt).kb = TissueDB(cnt).kf * TissueDB(cnt).M0b / (1 - TissueDB(cnt).M0b);
                
            case 'user'
                if isempty(user_tissue)
                    error('User tissue requested but UserTissue struct is empty.');
                end
                TissueDB(cnt) = user_tissue;
                if ~isfield(TissueDB(cnt),'Name'), TissueDB(cnt).Name = 'User'; end
                if ~isfield(TissueDB(cnt),'M0b'), error('User tissue must define M0b'); end
                if ~isfield(TissueDB(cnt),'T1f'), TissueDB(cnt).T1f = 1.0; end
                if ~isfield(TissueDB(cnt),'T2f'), TissueDB(cnt).T2f = 0.05; end
                if ~isfield(TissueDB(cnt),'T1b'), TissueDB(cnt).T1b = 1.0; end
                if ~isfield(TissueDB(cnt),'T2b'), TissueDB(cnt).T2b = 1e-5; end
                if ~isfield(TissueDB(cnt),'kf'), TissueDB(cnt).kf = 3.0; end
                if ~isfield(TissueDB(cnt),'kb')
                    TissueDB(cnt).kb = TissueDB(cnt).kf * TissueDB(cnt).M0b / (1 - TissueDB(cnt).M0b);
                end
            otherwise
                error('Unknown tissue name: %s', names{ii});
        end
    end
end

function TissueDB = verifyTissueParameters(TissueDB)
% Verify tissue parameter consistency and auto-correct if required
    for ti = 1:numel(TissueDB)
        params = TissueDB(ti);
        
        % Check M0b field existence
        if ~isfield(params,'M0b')
            error('MTCsim:MissingField', 'Tissue %d (%s) missing M0b field.', ti, params.Name);
        end
        
        % Check detailed balance: recompute desired kb
        kb_calc = params.kf * params.M0b / (1 - params.M0b);
        if ~isfield(params,'kb') || ~isnumeric(params.kb)
            TissueDB(ti).kb = kb_calc;
            warning('MTCsim:KBSet', 'kb not found or invalid for %s; set to %.6f.', params.Name, kb_calc);
        else
            tolerance = 1e-6;
            if abs(params.kb - kb_calc) > tolerance
                warning('MTCsim:DetailedBalance', ...
                    'Exchange rates for %s did not satisfy detailed balance. Overwriting kb (%.6g -> %.6g).', ...
                    params.Name, params.kb, kb_calc);
                TissueDB(ti).kb = kb_calc;
            end
        end
        
        % Check magnetization conservation
        if TissueDB(ti).M0b < 0 || TissueDB(ti).M0b > 1
            error('MTCsim:InvalidM0b', 'M0b for %s must be between 0 and 1', TissueDB(ti).Name);
        end
        
        % Check positive relaxation times
        reqFields = {'T1f','T2f','T1b','T2b','kf'};
        for f = 1:numel(reqFields)
            fld = reqFields{f};
            if ~isfield(TissueDB(ti), fld) || ~isnumeric(TissueDB(ti).(fld)) || ~(TissueDB(ti).(fld) > 0)
                error('MTCsim:InvalidRelaxation', 'Field %s must be positive for %s', fld, TissueDB(ti).Name);
            end
        end
    end
end

function [offsets, Noffsets, selected_offset, selected_idx] = setupOffsets(par)
% Setup offset frequencies for simulation
    if isempty(par.OffsetSweep)
        offsets = par.OffsetHz;  % single offset (scalar)
    else
        offsets = par.OffsetSweep(:).';
    end
    Noffsets = numel(offsets);
    
    % Determine selected offset
    if isempty(par.OffsetSweep)
        selected_offset = par.OffsetHz;
    else
        selected_offset = offsets(1);
    end
    
    % Find index of selected offset
    selected_idx = find(offsets == selected_offset, 1);
    if isempty(selected_idx)
        selected_idx = 1; 
        selected_offset = offsets(1);
    end
end

function [absorption, dispersion, linewidths_Hz] = computeLineshapes(TissueDB, offsets)
% Compute Lorentzian lineshapes for bound pool and their linewidths
    Ntissue = length(TissueDB);
    Noffsets = length(offsets);
    
    absorption = zeros(Ntissue, Noffsets);
    dispersion = zeros(Ntissue, Noffsets);
    linewidths_Hz = zeros(1, Ntissue);
    
    for ti = 1:Ntissue
        T2b = TissueDB(ti).T2b;
        
        % Calculate linewidth (FWHM of Lorentzian)
        linewidths_Hz(ti) = 1 / (pi * T2b); % Hz
        
        for io = 1:Noffsets
            Delta = 2*pi*offsets(io); % rad/s
            
            % Lorentzian lineshape components
            denominator = 1 + (Delta * T2b)^2;
            
            absorption(ti, io) = T2b / denominator;
            dispersion(ti, io) = (Delta * T2b^2) / denominator;
        end
        
        % Normalize absorption to maximum for comparison
        absorption(ti, :) = absorption(ti, :) / max(absorption(ti, :));
    end
end

function dM = blochMcConnellODE(M, omega1, Delta, params)
% Bloch-McConnell ODE for two-pool model - CORRECTED VERSION
% M: [Mx_f; My_f; Mz_f; Mx_b; My_b; Mz_b]
% omega1: RF amplitude (rad/s)
% Delta: off-resonance (rad/s)
% params: tissue parameters

    Mx_f = M(1); My_f = M(2); Mz_f = M(3);
    Mx_b = M(4); My_b = M(5); Mz_b = M(6);
    
    % Free pool Bloch equations in rotating frame (standard sign convention)
    dMx_f = -Delta * My_f - Mx_f/params.T2f;
    dMy_f =  Delta * Mx_f - My_f/params.T2f + omega1 * Mz_f;
    dMz_f = (1 - params.M0b - Mz_f)/params.T1f - omega1 * My_f ...
            - params.kf * Mz_f + params.kb * Mz_b;
    
    % Bound pool - assume instantaneous transverse saturation (short T2b)
    dMx_b = -Mx_b/params.T2b;
    dMy_b = -My_b/params.T2b;
    
    % CORRECTED: Bound pool saturation using Lorentzian absorption lineshape
    % Wb = π * ω1² * g(Δ) where g(Δ) is the absorption lineshape
    if omega1 ~= 0
        Wb = pi * (omega1^2) * params.T2b / (1 + (Delta * params.T2b)^2);
    else
        Wb = 0;
    end
    
    dMz_b = (params.M0b - Mz_b)/params.T1b - Wb * Mz_b ...
            + params.kf * Mz_f - params.kb * Mz_b;
    
    dM = [dMx_f; dMy_f; dMz_f; dMx_b; dMy_b; dMz_b];
end

function Mnext = rk4Step(Mcur, omega1, Delta, params, dt)
% 4th order Runge-Kutta integration step
    k1 = blochMcConnellODE(Mcur, omega1, Delta, params);
    k2 = blochMcConnellODE(Mcur + 0.5*dt*k1, omega1, Delta, params);
    k3 = blochMcConnellODE(Mcur + 0.5*dt*k2, omega1, Delta, params);
    k4 = blochMcConnellODE(Mcur + dt*k3, omega1, Delta, params);
    Mnext = Mcur + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end

function [M_final, Mz_free_curve, Mz_bound_curve, Mx_free_curve, My_free_curve, Mxy_free_curve] = ...
    integrateBlochMcConnell(M0, omega1_t, Delta, params, dt, tp, returnCurves)
% Integrate Bloch-McConnell equations through MT pulse
    if nargin < 7
        returnCurves = false;
    end
    
    N = length(omega1_t); % number of RF samples
    Mcur = M0;
    
    if returnCurves
        % store initial + after each RF sample => length N+1
        Mz_free_curve = zeros(N+1, 1);
        Mz_bound_curve = zeros(N+1, 1);
        Mx_free_curve = zeros(N+1, 1);
        My_free_curve = zeros(N+1, 1);
        Mxy_free_curve = zeros(N+1, 1);
        
        % Store initial conditions at index 1
        Mz_free_curve(1) = M0(3);
        Mz_bound_curve(1) = M0(6);
        Mx_free_curve(1) = M0(1);
        My_free_curve(1) = M0(2);
        Mxy_free_curve(1) = sqrt(M0(1)^2 + M0(2)^2);
    else
        Mz_free_curve = [];
        Mz_bound_curve = [];
        Mx_free_curve = [];
        My_free_curve = [];
        Mxy_free_curve = [];
    end
    
    % Use RK4 at each sample with local omega1; store result at n+1
    for n = 1:N
        omega1_now = omega1_t(n);
        Mcur = rk4Step(Mcur, omega1_now, Delta, params, dt);
        
        if returnCurves
            Mz_free_curve(n+1)  = Mcur(3);
            Mz_bound_curve(n+1) = Mcur(6);
            Mx_free_curve(n+1)  = Mcur(1);
            My_free_curve(n+1)  = Mcur(2);
            Mxy_free_curve(n+1) = sqrt(Mcur(1)^2 + Mcur(2)^2);
        end
    end
    
    M_final = Mcur;
end

function M_final = applyRecovery(M_start, Delta, params, dt, recoveryTime)
% Apply recovery period with no RF (omega1 = 0)
    steps = max(1, round(recoveryTime / dt));
    Mcur = M_start;
    for n = 1:steps
        Mcur = rk4Step(Mcur, 0, Delta, params, dt);
    end
    M_final = Mcur;
end

function MTR = computeMTR(M_final, params, alpha_deg)
% Compute Magnetization Transfer Ratio (simple single-shot readout model)
    alpha_rad = alpha_deg * pi/180;
    % Signal without MT pulse (equilibrium free pool longitudinal magnetization)
    S0 = sin(alpha_rad) * (1 - params.M0b);
    % Signal with MT pulse (use final free-pool longitudinal magnetization)
    S_MT = sin(alpha_rad) * M_final(3);
    if abs(S0) > eps
        MTR = (S0 - S_MT) / S0;
    else
        MTR = 0;
    end
end


%%% --- Plots
function createSummaryPlots(out, par, selected_offset)
% Create comprehensive summary plots with beautiful styling and optimized ranges
    if par.ShowLineshape
        fig = figure('Color','k','Position',[100 100 1400 1200]); % Larger for 6 subplots
        numSubplots = 6;
    else
        fig = figure('Color','k','Position',[100 100 1200 1000]);
        numSubplots = 5;
    end
    set(fig,'InvertHardcopy','off');

    t_ms = out.t*1e3;
    Nt = numel(out.TissueDB);
    colors = lines(Nt);
    subplotIdx = 1;
    
    % Enhanced color scheme with better contrast
    if Nt <= 3
        colors = [0.0, 0.8, 0.0;    % Green for Blood
                  0.8, 0.0, 0.0;    % Red for Myocardium  
                  0.9, 0.7, 0.0];   % Gold for Scar
    end
    
    % ==================== SUBPLOT 1: RF WAVEFORM ====================
    ax1 = subplot(numSubplots,1,subplotIdx); subplotIdx = subplotIdx + 1;
    B1_uT = out.B1_t * 1e6; % convert to microTesla
    plot(ax1, t_ms, B1_uT, 'y-', 'LineWidth',2.5, 'Color', [1, 0.8, 0.2]); % Golden yellow
    xlabel('Time (ms)','Color','w','FontSize',11,'FontWeight','bold');
    ylabel('B1 (\muT)','Color','w','FontSize',11,'FontWeight','bold');
    title(sprintf('Gaussian MT Pulse (Peak = %.1f \\muT, Flip_{eff} = %.1f°)', ...
          par.B1uT, out.flip_eff_deg),'Color','w','FontSize',12,'FontWeight','bold');
    set(ax1,'Color','k','XColor',[0.8,0.8,0.8],'YColor',[0.8,0.8,0.8],...
           'GridColor',[0.3,0.3,0.3],'GridAlpha',0.7);
    grid(ax1,'on');
    ylim(ax1, [0, par.B1uT * 1.2]); % Slightly above peak B1
    
    % ==================== SUBPLOT 2: Mz FREE ====================
    ax2 = subplot(numSubplots,1,subplotIdx); subplotIdx = subplotIdx + 1;
    hold(ax2,'on');
    for ti = 1:Nt
        plot(ax2, t_ms, out.Mz_free_curve(:,ti), 'LineWidth',2.2, 'Color', colors(ti,:));
    end
    xlabel('Time (ms)','Color','w','FontSize',11,'FontWeight','bold');
    ylabel('M_z (free)','Color','w','FontSize',11,'FontWeight','bold');
    title(sprintf('Free Pool Longitudinal Magnetization (Offset = %.0f Hz)', selected_offset),...
          'Color','w','FontSize',12,'FontWeight','bold');
    leg2 = legend(ax2, {out.TissueDB.Name}, 'Location','northeast',...
                 'TextColor','w','FontSize',10);
    set(leg2, 'Color', [0.1,0.1,0.1], 'TextColor', 'w', 'EdgeColor', [0.5,0.5,0.5]);
    set(ax2,'Color','k','XColor',[0.8,0.8,0.8],'YColor',[0.8,0.8,0.8],...
           'GridColor',[0.3,0.3,0.3],'GridAlpha',0.7);
    grid(ax2,'on');
    % BEAUTIFUL RANGE: Mz free typically from 0.6 to 1.0
    ylim(ax2, [0.7, 1.0]);
    
    % ==================== SUBPLOT 3: Mz BOUND ====================
    ax3 = subplot(numSubplots,1,subplotIdx); subplotIdx = subplotIdx + 1;
    hold(ax3,'on');
    for ti = 1:Nt
        plot(ax3, t_ms, out.Mz_bound_curve(:,ti),'--','LineWidth',2.0,'Color',colors(ti,:));
    end
    xlabel('Time (ms)','Color','w','FontSize',11,'FontWeight','bold');
    ylabel('M_z (bound)','Color','w','FontSize',11,'FontWeight','bold');
    title(sprintf('Bound Pool Longitudinal Magnetization (Offset = %.0f Hz)', selected_offset),...
          'Color','w','FontSize',12,'FontWeight','bold');
    leg3 = legend(ax3, {out.TissueDB.Name}, 'Location','northeast',...
                 'TextColor','w','FontSize',10);
    set(leg3, 'Color', [0.1,0.1,0.1], 'TextColor', 'w', 'EdgeColor', [0.5,0.5,0.5]);
    set(ax3,'Color','k','XColor',[0.8,0.8,0.8],'YColor',[0.8,0.8,0.8],...
           'GridColor',[0.3,0.3,0.3],'GridAlpha',0.7);
    grid(ax3,'on');
    % BEAUTIFUL RANGE: Mz bound typically from 0 to 0.25
    ylim(ax3, [0, 0.3]);
    
    % ==================== SUBPLOT 4: Mxy FREE ====================
    ax4 = subplot(numSubplots,1,subplotIdx); subplotIdx = subplotIdx + 1;
    hold(ax4,'on');
    for ti = 1:Nt
        plot(ax4, t_ms, out.Mxy_free_curve(:,ti), 'LineWidth',2.0, 'Color', colors(ti,:));
    end
    xlabel('Time (ms)','Color','w','FontSize',11,'FontWeight','bold');
    ylabel('M_{xy} (free)','Color','w','FontSize',11,'FontWeight','bold');
    title(sprintf('Free Pool Transverse Magnetization |M_{xy}| (Offset = %.0f Hz)', selected_offset),...
          'Color','w','FontSize',12,'FontWeight','bold');
    leg4 = legend(ax4, {out.TissueDB.Name}, 'Location','northeast',...
                 'TextColor','w','FontSize',10);
    set(leg4, 'Color', [0.1,0.1,0.1], 'TextColor', 'w', 'EdgeColor', [0.5,0.5,0.5]);
    set(ax4,'Color','k','XColor',[0.8,0.8,0.8],'YColor',[0.8,0.8,0.8],...
           'GridColor',[0.3,0.3,0.3],'GridAlpha',0.7);
    grid(ax4,'on');
    % BEAUTIFUL RANGE: Mxy typically very small, 0 to 0.02
   
    ylim(ax4, [0, 0.3]); % Standard range for small Mxy
    
    % ==================== SUBPLOT 5: BOUND POOL LINESHAPE ====================
    if par.ShowLineshape
        ax5 = subplot(numSubplots,1,subplotIdx); subplotIdx = subplotIdx + 1;
        hold(ax5,'on');
        
        % Plot absorption lineshape (solid) and dispersion (dashed)
        for ti = 1:Nt
            plot(ax5, out.offsets, out.lineshape_absorption(ti,:), '-', ...
                 'LineWidth',2.5, 'Color', colors(ti,:));
            plot(ax5, out.offsets, out.lineshape_dispersion(ti,:), '--', ...
                 'LineWidth',2.0, 'Color', colors(ti,:));
        end
        
        xlabel('Offset (Hz)','Color','w','FontSize',11,'FontWeight','bold');
        ylabel('Normalized Lineshape','Color','w','FontSize',11,'FontWeight','bold');
        title('Bound Pool Lorentzian Lineshape','Color','w','FontSize',12,'FontWeight','bold');
        subtitle('Absorption (solid) & Dispersion (dashed)','Color',[0.8,0.8,0.8],'FontSize',10);
        
        % Add text explaining the relationship
        text(0.02, 0.95, 'W_b(Δ) = πω₁² × g_{absorption}(Δ)', 'Units', 'normalized', ...
             'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', 'Parent', ax5);
        text(0.02, 0.88, 'Z-spectrum shape ∝ g_{absorption}(Δ)', 'Units', 'normalized', ...
             'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', 'Parent', ax5);
        
        leg5 = legend(ax5, {out.TissueDB.Name}, 'TextColor','w','Location','northeast',...
                     'FontSize',10);
        set(leg5, 'Color', [0.1,0.1,0.1], 'TextColor', 'w', 'EdgeColor', [0.5,0.5,0.5]);
        set(ax5,'Color','k','XColor',[0.8,0.8,0.8],'YColor',[0.8,0.8,0.8],...
               'GridColor',[0.3,0.3,0.3],'GridAlpha',0.7);
        grid(ax5,'on');
        % BEAUTIFUL RANGE: Lineshape from -0.5 to 1.1 to see both components
        ylim(ax5, [0. 0.4]);
    end

    % ==================== SUBPLOT 6: Z-SPECTRUM ====================
    ax6 = subplot(numSubplots,1,subplotIdx);
    hold(ax6,'on');
    
    % Create more distinctive markers and line styles
    markers = {'o','s','^','d','v','>','<'};
    lineStyles = {'-','--',':','-.'};
    
    for ti = 1:Nt
        marker = markers{mod(ti-1, length(markers)) + 1};
        lineStyle = lineStyles{mod(ti-1, length(lineStyles)) + 1};
        
        if par.PlotLogX
            semilogx(ax6, out.offsets, out.MTR_vs_offset(ti,:), ...
                    'LineStyle', lineStyle, 'Marker', marker, ...
                    'LineWidth',2.2, 'Color', colors(ti,:), ...
                    'MarkerSize',5, 'MarkerFaceColor', colors(ti,:));
        else
            plot(ax6, out.offsets, out.MTR_vs_offset(ti,:), ...
                 'LineStyle', lineStyle, 'Marker', marker, ...
                 'LineWidth',2.2, 'Color', colors(ti,:), ...
                 'MarkerSize',5, 'MarkerFaceColor', colors(ti,:));
        end
        
        % Overlay normalized lineshape to show correlation if enabled
        if par.ShowLineshape
            scaled_ls = out.lineshape_absorption(ti,:) * max(out.MTR_vs_offset(ti,:));
            if par.PlotLogX
                semilogx(ax6, out.offsets, scaled_ls, ':', ...
                        'LineWidth',1.8, 'Color', colors(ti,:));
            else
                plot(ax6, out.offsets, scaled_ls, ':', ...
                     'LineWidth',1.8, 'Color', colors(ti,:));
            end
        end
    end
    
    xlabel('Offset Frequency (Hz)','Color','w','FontSize',11,'FontWeight','bold');
    ylabel('MTR','Color','w','FontSize',11,'FontWeight','bold');
    
    if par.ShowLineshape
        title('Z-Spectrum vs Scaled Absorption Lineshape','Color','w','FontSize',12,'FontWeight','bold');
        subtitle('Solid: MTR, Dotted: Scaled Lineshape','Color',[0.8,0.8,0.8],'FontSize',10);
        % Create custom legend
        legendEntries = {};
        for ti = 1:Nt
            legendEntries{end+1} = [out.TissueDB(ti).Name ' MTR'];
            legendEntries{end+1} = [out.TissueDB(ti).Name ' lineshape'];
        end
        leg6 = legend(ax6, legendEntries, 'TextColor','w','Location','southeast',...
                     'NumColumns', 2, 'FontSize',9);
    else
        title('Magnetization Transfer Ratio (MTR) vs Offset','Color','w','FontSize',12,'FontWeight','bold');
        leg6 = legend(ax6, {out.TissueDB.Name}, 'TextColor','w','Location','southeast',...
                     'FontSize',10);
    end
    set(leg6, 'Color', [0.1,0.1,0.1], 'TextColor', 'w', 'EdgeColor', [0.5,0.5,0.5]);
    set(ax6,'Color','k','XColor',[0.8,0.8,0.8],'YColor',[0.8,0.8,0.8],...
           'GridColor',[0.3,0.3,0.3],'GridAlpha',0.7);
    grid(ax6,'on');
    % BEAUTIFUL RANGE: MTR typically from 0 to 0.6
    ylim(ax6, [0, 0.6]);

    % ==================== ENHANCE ALL PLOTS ====================
    % Link time-domain plots for synchronized zooming
    timeAxes = [ax1, ax2, ax3, ax4];
    linkaxes(timeAxes, 'x');
    
    % Set consistent time range for time-domain plots
    xlim(timeAxes, [min(t_ms), max(t_ms)]);
    
    % Add subtle background to legends for better readability
    allLegends = [leg2, leg3, leg4];
    if par.ShowLineshape, allLegends = [allLegends, leg5]; end
    allLegends = [allLegends, leg6];
    
    for leg = allLegends
        set(leg, 'Color', [0.05,0.05,0.05,0.9]); % Semi-transparent dark background
    end
    
    % Improve overall figure appearance
    set(fig, 'PaperPositionMode', 'auto');
    
    fprintf('Plot ranges optimized:\n');
    fprintf('  Mz_free: [0.6, 1.0]\n');
    fprintf('  Mz_bound: [0, 0.25]\n');
    fprintf('  Mxy_free: [0, 0.02]\n');
    fprintf('  MTR: [0, 0.6]\n');
end

function printSimulationSummary(out, TissueDB)
% Print comprehensive simulation summary
    fprintf('\n=== MTCsim Simulation Summary ===\n');
    fprintf('Pulse: %.2f μT Gaussian, %.2f ms, effective flip: %.3f deg\n', ...
        out.params.B1uT, out.params.PulseDur*1000, out.flip_eff_deg);
    fprintf('Offsets: %s Hz\n', mat2str(out.offsets));
    fprintf('Selected offset for detailed plots: %.0f Hz\n', out.selected_offset);
    
    if out.params.RecoveryTime > 0
        fprintf('Recovery time: %.2f ms after MT pulse\n', out.params.RecoveryTime*1000);
    end
    
    fprintf('\nTissue Results (MTR at first offset):\n');
    for ti = 1:length(TissueDB)
        fprintf('  %s: MTR = %.4f\n', TissueDB(ti).Name, out.MTR_vs_offset(ti,1));
    end
    
    fprintf('\nPeak transverse magnetization (Mxy) at selected offset:\n');
    for ti = 1:length(TissueDB)
        max_Mxy = max(out.Mxy_free_curve(:,ti));
        fprintf('  %s: max Mxy = %.6g\n', TissueDB(ti).Name, max_Mxy);
    end
    
    if out.params.ShowLineshape
        fprintf('\nBound Pool Lineshape Parameters:\n');
        for ti = 1:length(TissueDB)
            T2b = TissueDB(ti).T2b;
            linewidth_Hz = out.linewidths_Hz(ti);
            fprintf('  %s: T2b = %.2e s, Linewidth = %.0f Hz\n', ...
                TissueDB(ti).Name, T2b, linewidth_Hz);
        end
        fprintf('\nKey Physics: Z-spectrum shape follows bound pool absorption lineshape\n');
        fprintf('  W_b(Δ) = πω₁² × g_absorption(Δ)\n');
        fprintf('  MTR(Δ) ∝ g_absorption(Δ) × M0b × k_f\n');
    end
    fprintf('=== End Summary ===\n\n');
end



%{


% High-resolution simulation with lineshape analysis
out = MTCsim('OffsetSweep',0:100:6000, 'B1uT', 15, 'PulseSamples', 5000, ...
             'ShowLineshape', true, 'alpha', 30);

% Fast simulation without plots for batch processing  
out = MTCsim('OffsetSweep',[1000 2000 3000 4000 5000], 'ShowPlots', false);

% Custom tissue with specific lineshape properties
myTissue = struct('Name', 'Narrow Bound Pool', 'M0b', 0.12, ...
                  'T1f', 1.4, 'T2f', 0.06, 'T1b', 1.0, 'T2b', 2e-5, 'kf', 3.5);
out = MTCsim('Tissues', {'Myocardium', 'user'}, 'UserTissue', myTissue, ...
             'ShowLineshape', true);

% Study recovery effects
out = MTCsim('OffsetSweep', 500:250:5000, 'RecoveryTime', 0.1, 'B1uT', 12);



%}


