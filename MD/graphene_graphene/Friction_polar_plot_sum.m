% Plot directional-force polar plots for four \xi configurations

clear
close all

m = 22;
n = m + 1;
dirlist = {'1.0', '1.5', '2.0', '2.5', '4.0'};
xi_values = [1, 1.5, 2, 2.5, 4];

colors = [
    228, 26, 28;    % Red
    55, 126, 184;   % Blue
    77, 175, 74;    % Green
    152, 78, 163;   % Purple
    255, 127, 0     % Orange
    ] / 255;

for i = 1:length(dirlist)

    force = [];
    theta_deg_half = 0:3:180;
    ii = 1;

    for theta_deg = theta_deg_half
        filename = fullfile(dirlist{i}, sprintf('m%dn%d_tri30_%.1f_force_%d.txt', m, n, xi_values(i), theta_deg));

        if exist(filename, 'file')
            opts = detectImportOptions(filename);
            opts = setvartype(opts, {'double','double','double','double'});
            data = readtable(filename, opts);
            data = table2array(data);

            ForceX = data(:,3);
            ForceY = data(:,4);

            Forcei = ForceX.*cos(deg2rad(theta_deg)) + ForceY.*sin(deg2rad(theta_deg));
            force(ii) = max(abs(Forcei));
        else
            fprintf('Warning: File not found, skipping: %s\n', filename);
            force(ii) = NaN;
        end
        ii = ii + 1;
    end

    if ~isempty(force)
        force_full = [force, force(2:1:end-1)];
    else
        force_full = [];
    end
    alpha_deg_full = [theta_deg_half, 183:3:357];
    alpha_rad = deg2rad(alpha_deg_full);


    figure('Color', 'white','Position', [100, 100, 300, 400]); 

    polarplot(alpha_rad, force_full, 'o-', ...
        'LineWidth', 2.5, ...
        'MarkerSize', 6, ...
        'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:));
    ax = gca;
    ax.RAxis.Exponent = 0; 
    ax.RAxis.FontSize = 14;
    ax.RAxis.LineWidth = 1.2;

    ax.ThetaZeroLocation = 'right';
    ax.ThetaDir = 'counterclockwise';
    thetaticks(0:60:359);
    ax.FontSize=14;

    ax.GridColor = [0.8 0.8 0.8];
    ax.GridAlpha = 0.8;
    ax.FontName = 'Arial';

    q_val = 3 * xi_values(i)^2 / 4;

    title({['MD Simulation: $\xi = ', num2str(xi_values(i)), '$'], ...
        ['$q = ', num2str(q_val, '%.2f'), '$']}, ...
        'Interpreter', 'latex', 'FontSize', 20);

end