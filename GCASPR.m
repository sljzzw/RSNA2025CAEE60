function acquisition_order = GCASPR(Ny, Nz, M, num_interleaves, varargin)
    % SIMULATE_GCASPR Simulates the G-CASPR acquisition order in k-space.
    %
    % Inputs:
    %   Ny               - Number of phase encoding steps in y-direction (even).
    %   Nz               - Number of phase encoding steps in z-direction (even).
    %   M                - Number of rings (also lines per interleave).
    %   num_interleaves  - Number of spiral interleaves to simulate.
    %   plot_result      - (Optional) Boolean to control plotting. Default: true.
    %
    % Output:
    %   acquisition_order - (num_interleaves * M) x 2 matrix of [ky, kz] acquisition sequence.
    %
    % Assumptions:
    %   - Even Ny and Nz for centered k-space at (0,0).
    %   - Elliptical sampling region.
    %   - Golden angle ≈ 111.246 degrees.
    %   - Each interleave spirals outward with one full 360° turn.
    %
    % Author:
    %   Zhongwei Zhang, MD, PhD
    %   Washington University School of Medicine
    %   2025-10-11
    % 
    % Example:
    %   % Simulate for a 32x32 grid with 9 rings and 5 interleaves
    %   acquisition_order = GCASPR(32, 32, 9, 5);

    % Parse optional arguments
    p = inputParser;
    addOptional(p, 'plot_result', true, @islogical);
    parse(p, varargin{:});
    
    golden_angle = 180 * (sqrt(5) - 1);  % ≈111.246 degrees

    % Generate all grid points within the elliptical region (vectorized)
    ky_vals = -Ny/2 : Ny/2 - 1;
    kz_vals = -Nz/2 : Nz/2 - 1;
    [KY, KZ] = meshgrid(ky_vals, kz_vals);
    ellipse_mask = (KY / (Ny/2)).^2 + (KZ / (Nz/2)).^2 <= 1;
    points = [KY(ellipse_mask), KZ(ellipse_mask)];

    if isempty(points)
        error('No points in the elliptical region.');
    end

    % Compute normalized r for each point
    r = sqrt( (points(:,1) / (Ny/2)).^2 + (points(:,2) / (Nz/2)).^2 );

    % Divide into M rings with equal radial bins
    ring_bounds = linspace(0, 1, M+1);
    ring_points = cell(1, M);
    ring_theta = cell(1, M);
    
    for j = 1:M
        idx = (r >= ring_bounds(j)) & (r < ring_bounds(j+1));
        ring_points{j} = points(idx, :);
        if ~isempty(ring_points{j})
            thetas = atan2(ring_points{j}(:,2), ring_points{j}(:,1)) * 180 / pi;
            ring_theta{j} = mod(thetas, 360);  % 0 to 360
            % Sort by theta for reference
            [ring_theta{j}, sort_idx] = sort(ring_theta{j});
            ring_points{j} = ring_points{j}(sort_idx, :);
        else
            warning(['Ring ' num2str(j) ' has no points.']);
        end
    end

    % Simulate the acquisition order
    acquisition_order = zeros(num_interleaves * M, 2);
    acq_idx = 1;
    
    for i = 0:num_interleaves-1
        theta_start = mod(i * golden_angle, 360);
        for j = 1:M
            if isempty(ring_points{j})
                acquisition_order(acq_idx, :) = [0, 0];
            else
                theta_target = mod(theta_start + (j-1) / (M-1) * 360, 360);
                % Compute circular distance to theta_target
                diff1 = abs(ring_theta{j} - theta_target);
                diff2 = abs(ring_theta{j} - (theta_target + 360));
                diff3 = abs(ring_theta{j} - (theta_target - 360));
                circ_dist = min([diff1, diff2, diff3], [], 2);
                [~, min_idx] = min(circ_dist);
                selected_point = ring_points{j}(min_idx, :);
                acquisition_order(acq_idx, :) = selected_point;
            end
            acq_idx = acq_idx + 1;
        end
    end

    % Enhanced plotting with black background
    if p.Results.plot_result
        plot_gcaspr_results(acquisition_order, Ny, Nz, num_interleaves, M);
    end
end

function plot_gcaspr_results(acquisition_order, Ny, Nz, num_interleaves, M)
    % Create figure with black background
    fig = figure('Color', 'black', 'Position', [100, 100, 1200, 800]);
    
    % Main acquisition order plot
    subplot(2,3,[1,2,4,5]);
    hold on;
    
    % Create colormap for the spiral
    colors = parula(size(acquisition_order, 1));
    
    % Plot acquisition order with color progression and get handles for legend
    plot_handles = gobjects(1, 4); % Pre-allocate handles for legend
    
    % Plot the entire acquisition path as one entity for legend
    plot_handles(1) = plot(acquisition_order(:,1), acquisition_order(:,2), ...
        'Color', [0.2, 0.7, 1], 'LineWidth', 2, 'DisplayName', 'Acquisition Path');
    
    % Now plot with color progression (without legend entries)
    for i = 1:size(acquisition_order, 1)-1
        plot(acquisition_order(i:i+1, 1), acquisition_order(i:i+1, 2), ...
             'Color', colors(i,:), 'LineWidth', 2.5, ...
             'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', colors(i,:), ...
             'HandleVisibility', 'off');
    end
    
    % Highlight start and end points with proper legend handles
    plot_handles(2) = plot(acquisition_order(1,1), acquisition_order(1,2), 'go', ...
         'MarkerSize', 10, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'w', 'LineWidth', 2, ...
         'DisplayName', 'Start Point');
    
    plot_handles(3) = plot(acquisition_order(end,1), acquisition_order(end,2), 'ro', ...
         'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'w', 'LineWidth', 2, ...
         'DisplayName', 'End Point');
    
    % Plot elliptical boundary with proper handle
    theta = linspace(0, 2*pi, 100);
    ellipse_y = (Ny/2) * cos(theta);
    ellipse_z = (Nz/2) * sin(theta);
    plot_handles(4) = plot(ellipse_y, ellipse_z, 'w--', 'LineWidth', 1, ...
        'Color', [1,1,1,0.7], 'DisplayName', 'Elliptical Boundary');
    
    % Configure plot appearance
    axis equal;
    xlim([-Ny/2-1, Ny/2+1]);
    ylim([-Nz/2-1, Nz/2+1]);
    xlabel('k_y', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('k_z', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');
    title('G-CASPR Acquisition Order', 'Color', 'white', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Set axes properties for dark theme
    ax = gca;
    ax.Color = 'black';
    ax.XColor = 'white';
    ax.YColor = 'white';
    ax.GridColor = [0.5, 0.5, 0.5];
    ax.GridAlpha = 0.3;
    grid on;
    
    % Add colorbar to show acquisition progression
    c = colorbar;
    c.Color = 'white';
    c.Label.String = 'Acquisition Order';
    c.Label.Color = 'white';
    colormap(parula);
    
    % Create proper legend with correct handles
    legend(plot_handles, 'Location', 'best', ...
           'Color', [0, 0, 0], 'TextColor', 'white', 'EdgeColor', 'white', ...
           'FontSize', 10);
    
    % Time progression plot
    subplot(2,3,3);
    plot(1:size(acquisition_order,1), 'LineWidth', 3, 'Color', [0, 0.8, 1]);
    xlabel('Time Step', 'Color', 'white');
    ylabel('Point Index', 'Color', 'white');
    title('Temporal Progression', 'Color', 'white');
    grid on;
    set(gca, 'Color', 'black', 'XColor', 'white', 'YColor', 'white');
    
    % Radial distance progression
    subplot(2,3,6);
    radial_dist = sqrt(acquisition_order(:,1).^2 + acquisition_order(:,2).^2);
    plot(radial_dist, 'LineWidth', 3, 'Color', [0.8, 0.2, 0.8]);
    xlabel('Time Step', 'Color', 'white');
    ylabel('Radial Distance', 'Color', 'white');
    title('Radial Progression', 'Color', 'white');
    grid on;
    set(gca, 'Color', 'black', 'XColor', 'white', 'YColor', 'white');
    
    % Add information text box
    annotation('textbox', [0.02, 0.02, 0.4, 0.1], ...
               'String', sprintf('Ny=%d, Nz=%d\\nM=%d rings\\n%d interleaves\\nTotal points: %d', ...
               Ny, Nz, M, num_interleaves, size(acquisition_order,1)), ...
               'Color', 'white', 'BackgroundColor', [0, 0, 0], ...
               'FontSize', 10, 'Interpreter', 'none');
    
    hold off;
end