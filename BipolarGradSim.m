% BipolarGradSim.m
% DESCRIPTION:
%   Simulates bipolar gradient waveforms for Phase Contrast MRI (PC-MRI) velocity encoding.
%   The waveform consists of two trapezoidal lobes with opposite polarity, designed to null
%   the zeroth-order moment (M0) for stationary spins while encoding velocity via the
%   first-order moment (M1). The phase accumulation is computed using the equation:
%   ϕ(t) = ϕ₀ + γ∫G(τ)·r(τ)dτ, where r(τ) = r₀ + vτ + ½aτ².
%   Outputs include time, gradient waveform, moments (M0, M1, M2), phase accumulation,
%   and simulation parameters. Plots are generated to visualize results.
%
% SYNTAX:
%   [time, G, M0, M1, M2, phase, params] = BipolarGradSim()
%   [time, G, M0, M1, M2, phase, params] = BipolarGradSim('Parameter', value, ...)
%   [time, G, M0, M1, M2, phase, params] = BipolarGradSim('Gmax', 100)
%
% INPUT PARAMETERS:
%   'Gmax'        - Maximum gradient amplitude (mT/m), default: 50
%   'VENC'        - Velocity encoding (cm/s), default: 150 (adjusted for ±π phase)
%   'delta'       - Gradient lobe plateau duration (ms), default: 0.4
%   'T'           - Center-to-center lobe separation (ms), default: 0.4
%   'polarity'    - Polarity: -1 (negative-first), +1 (positive-first), default: -1
%   'rise_time'   - Gradient rise time (ms), default: 0.05
%   'fall_time'   - Gradient fall time (ms), default: 0.05
%   'dt'          - Time resolution (ms), default: 0.01
%   'velocities'  - Velocities to simulate (cm/s), default: [-120, -60, 0, 60, 120]
%   'r0'          - Initial position (m), default: 0
%   'accel'       - Acceleration (m/s²), default: 0
%   'phi0'        - Initial phase (rad), default: 0
%   'doPlot'      - Enable plotting, default: true
%   'style'       - Plot style: 'light' (default, white) or 'dark'
%
% OUTPUT PARAMETERS:
%   time       - Time vector (ms)
%   G          - Gradient waveform (T/m)
%   M0         - Zeroth moment: ∫G(τ)dτ (T·s/m)
%   M1         - First moment: ∫τG(τ)dτ (T·s²/m)
%   M2         - Second moment: ∫τ²G(τ)dτ (T·s³/m)
%   phase      - Phase accumulation for each velocity (rad)
%   params     - Structure containing simulation parameters and results
%
% GRADIENT MOMENTS:
%   M0: Position encoding moment (should be ~0 at end for stationary spins)
%   M1: Velocity encoding moment (primary for PC-MRI, ∝ v)
%   M2: Acceleration encoding moment (∝ a)
%
% PHASE EQUATION:
%   ϕ(t) = ϕ₀ + γ[r₀·M₀(t) + v·M₁(t) + ½a·M₂(t)]
%
% EXAMPLE USAGE:
%   % Default simulation
%   [t, G, M0, M1, M2, phase, params] = BipolarGradSim();
%
%   % Custom VENC and dark style
%   [time, G, M0, M1, M2, phase, params] = BipolarGradSim('Gmax', 30)
%
% REFERENCE:
%   Bernstein, M. A., King, K. F., & Zhou, X. J. (2004).
%   Handbook of MRI Pulse Sequences. Academic Press.
%
% AUTHOR: Zhongwei Zhang, MD, PhD
% Email: zhongweiz@wustl.edu
% Washington University School of Medicine


function [time, G, M0, M1, M2, phase, params] = BipolarGradSim(varargin)
%% Parse input parameters with validation
p = inputParser;
addParameter(p, 'Gmax', 50, @(x) isnumeric(x) && x > 0); % mT/m
addParameter(p, 'VENC', 150, @(x) isnumeric(x) && x > 0); % cm/s, adjusted for ±π
addParameter(p, 'delta', 0.4, @(x) isnumeric(x) && x > 0); % ms
addParameter(p, 'T', 0.4, @(x) isnumeric(x) && x > 0); % ms
addParameter(p, 'polarity', -1, @(x) isnumeric(x) && (x == -1 || x == 1));
addParameter(p, 'rise_time', 5e-3, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'fall_time', 5e-3, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'dt', 0.01, @(x) isnumeric(x) && x > 0);
addParameter(p, 'velocities', [-120, -60, 0, 60, 120], @isnumeric); % cm/s
addParameter(p, 'r0', 0, @isnumeric); % m
addParameter(p, 'accel', 0, @isnumeric); % m/s²
addParameter(p, 'phi0', 0, @isnumeric); % rad
addParameter(p, 'doPlot', true, @islogical);
addParameter(p, 'style', 'light', @(x) ischar(x) && ismember(x, {'light', 'dark'}));
parse(p, varargin{:});
params = p.Results;

%% Initialize output variables to ensure assignment
time = []; G = []; M0 = []; M1 = []; M2 = []; phase = []; params = p.Results;

%% Define physical constants
gamma = 2.675e8; % Gyromagnetic ratio for 1H (rad/T·s), precise value

%% Convert units to SI for calculations
Gmax_T = params.Gmax * 1e-3; % mT/m to T/m
delta_s = params.delta * 1e-3; % ms to s
T_s = params.T * 1e-3; % ms to s
rise_time_s = params.rise_time * 1e-3; % ms to s
fall_time_s = params.fall_time * 1e-3; % ms to s
dt_s = params.dt * 1e-3; % ms to s
VENC_mps = params.VENC * 1e-2; % cm/s to m/s
velocities_mps = params.velocities * 1e-2; % cm/s to m/s
r0_m = params.r0; % m
accel_mps2 = params.accel; % m/s²
phi0_rad = params.phi0; % rad

%% Validate timing parameters
if dt_s <= 0
    error('Time resolution (dt) must be positive.');
end
if delta_s <= 0
    error('Plateau duration (delta) must be positive.');
end
if T_s <= 0
    error('Center-to-center separation (T) must be positive.');
end
if rise_time_s + fall_time_s + delta_s <= 0
    error('Single lobe duration (rise_time + delta + fall_time) must be positive.');
end

%% Calculate gap between lobes
% Gap is the time between the end of the first lobe and the start of the second
gap_s = T_s - delta_s; % T is center-to-center, so gap adjusts for plateau duration
if gap_s < 0
    warning('T < delta; setting gap to 0 to avoid negative duration.');
    gap_s = 0;
end

%% Create time vector
% Total duration includes both lobes and gap
total_time = rise_time_s + delta_s + fall_time_s + gap_s + rise_time_s + delta_s + fall_time_s;
if total_time <= 0
    error('Total waveform duration is non-positive. Check timing parameters.');
end
time_s = 0:dt_s:total_time; % Time vector in seconds
if isempty(time_s)
    error('Time vector is empty. Check dt and total_time.');
end
time_ms = time_s * 1e3; % Convert to ms for output


%% Generate bipolar gradient waveform
G = zeros(size(time_s)); % Initialize gradient waveform (T/m)

% Define time segments for the trapezoidal lobes
t1 = rise_time_s; % End of first rise
t2 = t1 + delta_s; % End of first plateau
t3 = t2 + fall_time_s + gap_s; % Start of second rise
t4 = t3 + rise_time_s; % End of second rise
t5 = t4 + delta_s; % End of second plateau
t6 = t5 + fall_time_s; % End of second fall

% First lobe (rise, plateau, fall)
% Rise: Linear ramp from 0 to ±Gmax
idx1 = time_s <= t1;
if rise_time_s > 0
    G(idx1) = params.polarity * Gmax_T * (time_s(idx1) / rise_time_s);
end
% Plateau: Constant at ±Gmax
idx2 = (time_s > t1) & (time_s <= t2);
G(idx2) = params.polarity * Gmax_T;
% Fall: Linear ramp from ±Gmax to 0
idx3 = (time_s > t2) & (time_s <= t2 + fall_time_s);
if fall_time_s > 0
    fall_segment = time_s(idx3) - t2;
    G(idx3) = params.polarity * Gmax_T * (1 - fall_segment / fall_time_s);
end

% Gap: G = 0 between lobes

% Second lobe (opposite polarity)
% Rise: Linear ramp from 0 to ∓Gmax
idx4 = (time_s > t3) & (time_s <= t4);
if rise_time_s > 0
    rise_segment = time_s(idx4) - t3;
    G(idx4) = -params.polarity * Gmax_T * (rise_segment / rise_time_s);
end
% Plateau: Constant at ∓Gmax
idx5 = (time_s > t4) & (time_s <= t5);
G(idx5) = -params.polarity * Gmax_T;
% Fall: Linear ramp from ∓Gmax to 0
idx6 = (time_s > t5) & (time_s <= t6);
if fall_time_s > 0
    fall_segment = time_s(idx6) - t5;
    G(idx6) = -params.polarity * Gmax_T * (1 - fall_segment / fall_time_s);
end

% Trim vectors to exact waveform duration
valid_idx = time_s <= t6;
if ~any(valid_idx)
    error('No valid time points after trimming. Check timing parameters.');
end
time_s = time_s(valid_idx);
time_ms = time_ms(valid_idx);
G = G(valid_idx);
time = time_ms;

%% Calculate gradient moments using trapezoidal integration
% M0 = ∫ G(τ) dτ (T·s/m)
M0 = zeros(size(time_s));
for i = 2:length(time_s)
    f_prev = G(i-1);
    f_curr = G(i);
    M0(i) = M0(i-1) + ((f_prev + f_curr) / 2) * dt_s;
end

% M1 = ∫ τ G(τ) dτ (T·s²/m)
M1 = zeros(size(time_s));
for i = 2:length(time_s)
    f_prev = time_s(i-1) * G(i-1);
    f_curr = time_s(i) * G(i);
    M1(i) = M1(i-1) + ((f_prev + f_curr) / 2) * dt_s;
end

% M2 = ∫ τ² G(τ) dτ (T·s³/m)
M2 = zeros(size(time_s));
for i = 2:length(time_s)
    f_prev = time_s(i-1)^2 * G(i-1);
    f_curr = time_s(i)^2 * G(i);
    M2(i) = M2(i-1) + ((f_prev + f_curr) / 2) * dt_s;
end

%% Calculate phase accumulation
% Phase = ϕ₀ + γ [r₀·M₀(t) + v·M₁(t) + ½a·M₂(t)]
phase = zeros(length(velocities_mps), length(time_s));
for v = 1:length(velocities_mps)
    phase(v,:) = phi0_rad + gamma * (r0_m * M0 + velocities_mps(v) * M1 + 0.5 * accel_mps2 * M2);
end

%% Store parameters and expected M1
% Expected M1 for phase = π at VENC: M1 = π / (γ * VENC)
params.gamma = gamma;
params.VENC_mps = VENC_mps;
params.velocities_mps = velocities_mps;
params.expected_M1 = pi / (gamma * VENC_mps); % Corrected for ±π phase at VENC

%% Plot results if requested
if params.doPlot
    plotResults(time_ms, G*1e3, M0, M1, M2, phase, params);
end

end

%% Plotting function
function plotResults(time_ms, G_mTm, M0, M1, M2, phase, params)
% Creates a 5-subplot figure showing the gradient waveform, moments, and phase
% accumulation, with a parameter summary on the right.

    % Initialize figure
    figure('Position', [100, 100, 1600, 900], 'Color', [1, 1, 1]);
    
    % Set global font properties
    set(findall(gcf,'-property','FontSize'),'FontSize',16);
    set(findall(gcf,'-property','FontName'),'FontName','Arial');
    set(findall(gcf,'-property','LineWidth'),'LineWidth',2);
    
    % Set plot style
    if strcmp(params.style, 'dark')
        bg_color = [0.15, 0.15, 0.15]; % Dark gray
        text_color = [1, 1, 1]; % White text
        grid_color = [0.6, 0.6, 0.6]; % Light gray grid
    else
        bg_color = [1, 1, 1]; % White
        text_color = [0, 0, 0]; % Black text
        grid_color = [0.8, 0.8, 0.8]; % Light gray grid
    end
    
    % Convert phase to degrees
    phase_deg = rad2deg(phase);
    
    % Define subplot height
    subplot_height = 0.12;
    
    % Subplot 1: Bipolar Gradient Waveform
    subplot(5,1,1, 'Position', [0.1, 0.82, 0.65, subplot_height]);
    plot(time_ms, G_mTm, 'b-', 'LineWidth', 2);
    hold on;
    
    % Mark key time points
    lobe_duration = params.delta + params.rise_time + params.fall_time;
    key_times = [0, params.rise_time, params.rise_time + params.delta, ...
                 params.rise_time + params.delta + params.fall_time, ...
                 lobe_duration + params.T - params.delta, ...
                 lobe_duration + params.T - params.delta + params.rise_time, ...
                 lobe_duration + params.T - params.delta + params.rise_time + params.delta, ...
                 lobe_duration + params.T - params.delta + params.rise_time + params.delta + params.fall_time];
    key_times = unique(key_times(key_times <= max(time_ms)));
    plot(key_times, zeros(size(key_times)), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    
    set(gca, 'Color', bg_color, 'XColor', text_color, 'YColor', text_color);
    grid on; set(gca, 'GridColor', grid_color, 'GridAlpha', 0.3);
    xlabel('Time (ms)', 'Color', text_color);
    ylabel('Gradient (mT/m)', 'Color', text_color);
    % Set title based on polarity
    if params.polarity == -1
        title_str = 'Bipolar Gradient: Negative-first [-G, +G]';
    else
        title_str = 'Bipolar Gradient: Positive-first [+G, -G]';
    end
    title(title_str, 'Color', text_color);
    legend('Gradient', 'Key points', 'Location', 'northwest', 'TextColor', text_color, ...
           'Color', [1, 1, 1], 'EdgeColor', 'none', 'FontSize', 9);
    
    % Subplot 2: Zeroth-Order Moment (M0)
    subplot(5,1,2, 'Position', [0.1, 0.63, 0.65, subplot_height]);
    plot(time_ms, M0, 'm-', 'LineWidth', 2);
    hold on;
    plot([min(time_ms), max(time_ms)], [0, 0], 'r--', 'LineWidth', 1);
    set(gca, 'Color', bg_color, 'XColor', text_color, 'YColor', text_color);
    grid on; set(gca, 'GridColor', grid_color, 'GridAlpha', 0.3);
    xlabel('Time (ms)', 'Color', text_color);
    ylabel('M0 (T·s/m)', 'Color', text_color);
    title('Zeroth-Order Moment', 'Color', text_color);
    legend('M0(t)', 'Ideal (0)', 'Location', 'northwest', 'TextColor', text_color, ...
           'Color', [1, 1, 1], 'EdgeColor', 'none', 'FontSize', 9);
    
    % Subplot 3: First-Order Moment (M1)
    subplot(5,1,3, 'Position', [0.1, 0.44, 0.65, subplot_height]);
    plot(time_ms, M1, 'k-', 'LineWidth', 2);
    hold on;
    expected_M1 = params.expected_M1 * sign(M1(end));
    plot([min(time_ms), max(time_ms)], [expected_M1, expected_M1], 'r--', 'LineWidth', 1);
    set(gca, 'Color', bg_color, 'XColor', text_color, 'YColor', text_color);
    grid on; set(gca, 'GridColor', grid_color, 'GridAlpha', 0.3);
    xlabel('Time (ms)', 'Color', text_color);
    ylabel('M1 (T·s²/m)', 'Color', text_color);
    title('First-Order Moment', 'Color', text_color);
    legend('M1(t)', 'Theoretical M1', 'Location', 'northwest', 'TextColor', text_color, ...
           'Color', [1, 1, 1], 'EdgeColor', 'none', 'FontSize', 9);
    
    % Subplot 4: Second-Order Moment (M2)
    subplot(5,1,4, 'Position', [0.1, 0.25, 0.65, subplot_height]);
    plot(time_ms, M2, 'c-', 'LineWidth', 2);
    set(gca, 'Color', bg_color, 'XColor', text_color, 'YColor', text_color);
    grid on; set(gca, 'GridColor', grid_color, 'GridAlpha', 0.3);
    xlabel('Time (ms)', 'Color', text_color);
    ylabel('M2 (T·s³/m)', 'Color', text_color);
    title('Second-Order Moment', 'Color', text_color);
    legend('M2(t)', 'Location', 'northwest', 'TextColor', text_color, ...
           'Color', [1, 1, 1], 'EdgeColor', 'none', 'FontSize', 9);
    
    % Subplot 5: Phase Accumulation
    subplot(5,1,5, 'Position', [0.1, 0.06, 0.65, subplot_height]);
    colors = lines(length(params.velocities));
    for v = 1:length(params.velocities)
        plot(time_ms, phase_deg(v,:), 'Color', colors(v,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('v = %d cm/s', params.velocities(v)));
        hold on;
    end
    % Add reference lines for ±π (180°) at VENC
    plot([min(time_ms), max(time_ms)], [180, 180], 'r--', 'LineWidth', 1, 'DisplayName', 'VENC (+π)');
    plot([min(time_ms), max(time_ms)], [-180, -180], 'r--', 'LineWidth', 1, 'DisplayName', 'VENC (-π)');
    set(gca, 'Color', bg_color, 'XColor', text_color, 'YColor', text_color);
    grid on; set(gca, 'GridColor', grid_color, 'GridAlpha', 0.3);
    xlabel('Time (ms)', 'Color', text_color);
    ylabel('Phase (deg)', 'Color', text_color);
    title('Phase Accumulation', 'Color', text_color);
    legend('Location', 'northwest', 'TextColor', text_color, 'Color', [1, 1, 1], ...
           'EdgeColor', 'none', 'FontSize', 9, 'NumColumns', 2);
    
    % Parameter annotation box
    annotation('textbox', [0.77, 0.06, 0.21, 0.4], 'String', ...
        sprintf(['Parameters:\n' ...
                'G_{max} = %.1f mT/m\n' ...
                'VENC = %.1f cm/s\n' ...
                'δ = %.1f ms\n' ...
                'T = %.1f ms\n' ...
                'Polarity = %d\n' ...
                'r_0 = %.2e m\n' ...
                'Accel = %.2e m/s²\n' ...
                'φ_0 = %.2f rad\n' ...
                'M0_{end} = %.2e T·s/m\n' ...
                'M1_{end} = %.2e T·s²/m\n' ...
                'M1_{exp} = %.2e T·s²/m\n' ...
                'M2_{end} = %.2e T·s³/m'], ...
                params.Gmax, params.VENC, params.delta, params.T, params.polarity, ...
                params.r0, params.accel, params.phi0, M0(end), M1(end), abs(params.expected_M1), M2(end)), ...
        'BackgroundColor', [1, 1, 1], 'EdgeColor', 'black', 'FontSize', 10, ...
        'FitBoxToText', 'on', 'Margin', 3);
end