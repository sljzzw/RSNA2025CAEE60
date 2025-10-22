function sim_out = T2Prep_MLEV4(varargin)
% T2Prep_MLEV4 - Simulate T2 preparation with MLEV-4 refocusing scheme
%
% Syntax:
%   sim_out = T2Prep_MLEV4('T1', [1000 1500], 'T2', [50 250], 'TE', 50, ...
%                          'dt', 0.1, 'M0', 1, 'spoiler', true, ...
%                          'spoiler_grad', 10, 'doPlot', true)
%
% Description:
%   Simulate a simple Bloch-style relaxation model of a T2-preparation
%   module using an MLEV-4 refocusing pattern (+X, +Y, +Y, +X).
%   RF pulses are applied instantaneously at their scheduled times;
%   relaxation (T1, T2) is applied between events with time step dt.
%
% Inputs (name-value):
%   'T1'         - vector of T1 values in ms (default 1200)
%   'T2'         - vector of T2 values in ms (default 45)
%   'TE'         - total T2-prep duration in ms (default 80)
%   'dt'         - simulation time-step in ms (default 0.1)
%   'M0'         - equilibrium magnetization (default 1)
%   'spoiler'    - logical, whether to apply spoiling after restore (default true)
%   'spoiler_grad'- spoiler gradient strength (not used in simple model) (default 10)
%   'doPlot'     - logical, whether to produce figures (default true)
%
% Output:
%   sim_out - struct with fields:
%             time, Mz, Mxy, Mx, My, RFtimes, RFindices, Mz_end, Mxy_end
%
% Author:
%   Zhongwei Zhang, MD, PhD
%   Washington University School of Medicine
%   2025-10-11
%
%{

% The simulation used tissue parameters for myocardium and arterial blood

T1 = [950, 1500]; 
T2 = [55, 250]; 
TE = 60; 
dt = 0.1; 
M0 = 1; 
spoiler = true; 
spoiler_grad = 10;

sim_out = T2Prep_MLEV4('T1', T1, 'T2', T2, 'TE', TE, 'dt', dt, ...
                       'M0', M0, 'spoiler', spoiler, ...
                       'spoiler_grad', spoiler_grad, ...
                       'doPlot', true);

fprintf('Tissue 1 - Final Longitudinal Magnetization (Mz_end): %.4f\n', sim_out.Mz_end(1));
fprintf('Tissue 1 - Final Transverse Magnetization (Mxy_end): %.4f\n', sim_out.Mxy_end(1));
fprintf('Tissue 2 - Final Longitudinal Magnetization (Mz_end): %.4f\n', sim_out.Mz_end(2));
fprintf('Tissue 2 - Final Transverse Magnetization (Mxy_end): %.4f\n', sim_out.Mxy_end(2));


%}



% ---------------------------
% Input parsing & validation
% ---------------------------
p = inputParser;
addParameter(p,'T1',1200,@(x)isnumeric(x) && all(x>0));
addParameter(p,'T2',45,@(x)isnumeric(x) && all(x>0));
addParameter(p,'TE',80,@(x)isnumeric(x) && isscalar(x) && (x>0));
addParameter(p,'dt',0.1,@(x)isnumeric(x) && isscalar(x) && (x>0));
addParameter(p,'M0',1,@isnumeric);
addParameter(p,'spoiler',true,@islogical);
addParameter(p,'spoiler_grad',10,@isnumeric);
addParameter(p,'doPlot',true,@islogical);
parse(p,varargin{:});

T1 = p.Results.T1(:);
T2 = p.Results.T2(:);
TE = p.Results.TE;
dt = p.Results.dt;
M0 = p.Results.M0;
spoiler = p.Results.spoiler;
spoiler_grad = p.Results.spoiler_grad;
doPlot = p.Results.doPlot;

if numel(T1) ~= numel(T2)
    error('T1 and T2 must be the same length.');
end
nTissue = numel(T1);

% ---------------------------
% RF scheme (MLEV-4)
% ---------------------------
nRef = 4;
rf_exc = 90;                          % excitation
rf_ref = repmat(180,1,nRef);          % 4 refocusing pulses
rf_restore = 90;                      % restore
rf_phase_exc = 0;
rf_phase_ref = [0 90 90 0];           % MLEV-4: +X, +Y, +Y, +X
rf_phase_restore = 180;

% timings
t_exc = 0;
t_restore = TE;
ref_interval = TE / nRef;
t_ref = (t_exc + ref_interval/2) : ref_interval : (t_restore - ref_interval/2);
rf_events = [t_exc, t_ref, t_restore];
flip_angles = [rf_exc, rf_ref, rf_restore];
flip_phases = [rf_phase_exc, rf_phase_ref, rf_phase_restore];

% ---------------------------
% Time vector that includes RF events exactly
% ---------------------------
time_base = 0:dt:(TE + 10);                      % simulate a bit after restore
time = unique(sort([time_base, rf_events]));
Nt = numel(time);

% map event times to indices in time vector
rf_indices = zeros(1,numel(rf_events));
for k = 1:numel(rf_events)
    [~, rf_indices(k)] = min(abs(time - rf_events(k)));
end

% ---------------------------
% Initialize magnetization arrays
% ---------------------------
Mz = zeros(nTissue, Nt); Mz(:,1) = M0;
Mx = zeros(nTissue, Nt); My = zeros(nTissue, Nt);
Mxy = zeros(nTissue, Nt);
Mxy(:,1) = sqrt(Mx(:,1).^2 + My(:,1).^2);

% ---------------------------
% Main simulation loop
% ---------------------------
for ti = 1:Nt-1
    t_current = time(ti);
    t_next = time(ti+1);
    dt_step = t_next - t_current;

    % current vectors
    Mz_curr = Mz(:,ti);
    Mx_curr = Mx(:,ti);
    My_curr = My(:,ti);

    % apply instantaneous RF if this index is an RF event
    if ismember(ti, rf_indices)
        rf_idx = find(rf_indices == ti, 1);
        for tt = 1:nTissue
            [Mx_curr(tt), My_curr(tt), Mz_curr(tt)] = ...
                applyRF(Mx_curr(tt), My_curr(tt), Mz_curr(tt), ...
                        flip_angles(rf_idx), flip_phases(rf_idx));
        end
    end

    % relaxation over dt_step
    E1 = exp(-dt_step ./ T1);
    E2 = exp(-dt_step ./ T2);

    Mz(:,ti+1) = Mz_curr .* E1 + M0 .* (1 - E1);
    Mx(:,ti+1) = Mx_curr .* E2;
    My(:,ti+1) = My_curr .* E2;
    Mxy(:,ti+1) = sqrt(Mx(:,ti+1).^2 + My(:,ti+1).^2);
end

% ---------------------------
% Optional spoiling after restore
% ---------------------------
if spoiler
    [~, idx_te] = min(abs(time - TE));
    if idx_te <= Nt
        Mx(:,idx_te:end) = 0;
        My(:,idx_te:end) = 0;
        Mxy(:,idx_te:end) = 0;
    end
end

% ---------------------------
% Prepare outputs
% ---------------------------
sim_out.time = time;
sim_out.Mz = Mz;
sim_out.Mxy = Mxy;
sim_out.Mx = Mx;
sim_out.My = My;
sim_out.RFtimes = rf_events;
sim_out.RFindices = rf_indices;
sim_out.Mz_end = Mz(:,end);
sim_out.Mxy_end = Mxy(:,end);

%% Plotting
if doPlot
    fig = figure('Color','k','Position',[100 100 1000 700],'InvertHardcopy','off');
    colors = lines(nTissue);

    % ---------- Subplot 1 ----------
    ax1 = subplot(3,1,1); hold(ax1,'on');
    c_ref = [1 0 0; 0 1 1; 1 0 0; 0 1 1];
    for k = 1:numel(rf_events)
        if k==1
            c = [0 1 0]; lbl='Excitation';
        elseif k==numel(rf_events)
            c = [1 0 1]; lbl='Restore';
        else
            c = c_ref(k-1,:); lbl='Refocusing';
        end
        stem(ax1, rf_events(k), flip_angles(k), 'filled', ...
            'LineWidth', 2, 'MarkerSize', 7, 'Color', c);
        text(ax1, rf_events(k), flip_angles(k)+18, ...
            sprintf('%d° (%s)', flip_angles(k), lbl), ...
            'HorizontalAlignment','center','FontWeight','bold', ...
            'Color','w','FontSize',9);
    end
    % Hide y axis
    ax1.YColor = 'none';
    ax1.YTick = [];
    ylabel(ax1,'','Color','w');
    title(ax1,'MLEV-4 RF Pulse Sequence','Color','w','FontWeight','bold');
    grid(ax1,'on');
    xlim(ax1,[0 TE+10]);
    set(ax1,'XTick',rf_events,'XTickLabel',rf_events, ...
        'Color','k','XColor','w','GridColor',[0.4 0.4 0.4]);

    % ---------- Subplot 2 ----------
    ax2 = subplot(3,1,2); hold(ax2,'on');
    for tt = 1:nTissue
        plot(ax2, time, Mz(tt,:), 'LineWidth', 2, 'Color', colors(tt,:));
    end
    for k = 1:numel(rf_events)
        xline(ax2, rf_events(k),'--','LineWidth',1,'Alpha',0.7,'Color',[0.5 0.5 0.5]);
    end
    ylabel(ax2,'M_z','Color','w');
    title(ax2,'Longitudinal Magnetization (Mz)','Color','w','FontWeight','bold');
    grid(ax2,'on');
    xlim(ax2,[0 TE+10]);
    ylim(ax2,[-0.2 1]);
    set(ax2,'XTick',rf_events,'XTickLabel',rf_events, ...
        'Color','k','XColor','w','YColor','w','GridColor',[0.4 0.4 0.4]);
    legend(arrayfun(@(i) sprintf('T1=%.0f ms, T2=%.0f ms', T1(i), T2(i)), ...
        1:nTissue,'UniformOutput',false), ...
        'TextColor','w','Color','none','EdgeColor','none','Location','best');

    % ---------- Subplot 3 ----------
    ax3 = subplot(3,1,3); hold(ax3,'on');
    for tt = 1:nTissue
        plot(ax3, time, Mxy(tt,:), 'LineWidth', 2, 'Color', colors(tt,:));
    end
    for k = 1:numel(rf_events)
        xline(ax3, rf_events(k),'--','LineWidth',1,'Alpha',0.7,'Color',[0.5 0.5 0.5]);
    end
    xlabel(ax3,'Time (ms)','Color','w');
    ylabel(ax3,'|M_{xy}|','Color','w');
    title(ax3,'Transverse Magnetization (Mxy)','Color','w','FontWeight','bold');
    grid(ax3,'on');
    xlim(ax3,[0 TE+10]);
    set(ax3,'XTick',rf_events,'XTickLabel',rf_events, ...
        'Color','k','XColor','w','YColor','w','GridColor',[0.4 0.4 0.4]);
    legend(arrayfun(@(i) sprintf('T1=%.0f ms, T2=%.0f ms', T1(i), T2(i)), ...
        1:nTissue,'UniformOutput',false), ...
        'TextColor','w','Color','none','EdgeColor','none','Location','best');

    sgtitle(sprintf('T2 Preparation with MLEV-4 (TE = %d ms)', TE), ...
        'Color','w','FontSize',13,'FontWeight','bold');
end

end

%% RF helper function
function [Mx_out,My_out,Mz_out] = applyRF(Mx,My,Mz,flip_deg,phase_deg)
% applyRF - rotate the magnetization vector by flip angle about an axis
% defined by phase_deg in the transverse (xy) plane.
%
% This uses the rotation: R = Rz(phase) * Ry(flip) * Rz(-phase).
flip_rad = deg2rad(flip_deg);
phase_rad = deg2rad(phase_deg);

cp = cos(phase_rad);
sp = sin(phase_rad);
cf = cos(flip_rad);
sf = sin(flip_rad);

% Compose rotation matrix explicitly (3x3)
R = [ cp^2*(1-cf)+cf,      cp*sp*(1-cf),    -sp*sf;
      cp*sp*(1-cf),   sp^2*(1-cf)+cf,     cp*sf;
      sp*sf,          -cp*sf,             cf ];

M_vec = [Mx; My; Mz];
M_rot = R * M_vec;

Mx_out = M_rot(1);
My_out = M_rot(2);
Mz_out = M_rot(3);
end
