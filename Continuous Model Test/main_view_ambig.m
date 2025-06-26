%% Generate figures
fig_num = 1;
font_val = 16;
line_val = 2;

% clear; clc; close all;
section = 1;
% Demonstrator settings
shapes = ["rrc"];
shape_names = ["RRC"];
truncation = false;
res = 100;
alpha = 1;
q = 2;

% Set time and frequency ranges
T = 1;
delta_f = 1/T;
t_range = linspace(-3*T,3*T,res);
f_range = [0,1,2] * delta_f;

% Convert to meshes
if not(truncation)
    t_range_mat = t_range.' .* ones(length(t_range),length(f_range));
    f_range_mat = f_range .* ones(length(t_range),length(f_range));
end

% Create matrix of results
result_2d = zeros(length(t_range),length(f_range),length(shapes));
for k = 1:length(shapes)
    if truncation
        for m = 1:length(t_range)
            for n = 1:length(f_range)
                result_2d(m,n,k) = ambig_direct(t_range(m),f_range(n),T,shapes(k),alpha,q,res);
            end
        end
    else
        result_2d(:,:,k) = ambig(t_range_mat,f_range_mat,T,shapes(k),alpha);
    end
end

legend_vec = ["r-","k--","b-."];

% Plot a figure for each shape
for k = 1:length(shapes)
    figure(k+length(shapes)*(section-1));
    hold on;
    for n = 1:length(f_range)
        plot(t_range,abs(result_2d(:,n,k)),legend_vec(n));
    end
    grid on;
    % title(sprintf("Amplitude of ambiguity function using %s pulses",shape_names(k)))
    xlabel("t (in T_s)")
    ylabel(sprintf("|A_{%s}(t,f)|",shapes(k)))
    legend("f=0\Deltaf","f=1\Deltaf","f=2\Deltaf")

    % Set font size
    set(gca, 'FontSize', font_val);

    % Adjust plotted line and border line size
    lines = findall(gcf, 'Type', 'Line');
    set(lines, 'LineWidth', line_val);
    ax = gca;
    ax.LineWidth = line_val;

    % Save figure
    saveas(figure(fig_num), sprintf('Figure%d_%s.eps',fig_num,string(datetime('now', 'format', 'MM.dd.uuuu_HH.mm'))), 'epsc');
end

%% Generate 3D view
section = 2;
% Set time and frequency ranges
t_range = linspace(-4*T,4*T,res);
f_range = linspace(-4/T,4/T,res);

% Convert to meshes
t_range_mat = t_range.' .* ones(res,length(f_range));
f_range_mat = f_range .* ones(res,length(f_range));

% Create matrix of results
result_3d = zeros(length(t_range),length(f_range),length(shapes));
for k = 1:length(shapes)
    if truncation
        for m = 1:length(t_range)
            for n = 1:length(f_range)
                result_3d(m,n,k) = ambig_direct(t_range(m),f_range(n),T,shapes(k),alpha,q,res);
            end
        end
    else
        result_3d(:,:,k) = ambig(t_range_mat,f_range_mat,T,shapes(k),alpha);
    end
end

% Plot a 3D mesh for each shape
for k = 1:length(shapes)
    figure(k+length(shapes)*(section-1));
    mesh(t_range,f_range,abs(result_3d(:,:,k)));
    title(sprintf("Amplitude of ambiguity function using %s pulses",shapes(k)))
    ylabel("t (in T)")
    xlabel("f (in 1/T)")
    zlabel(sprintf("|A_{%s}(t,f)|",shapes(k)))
end
