function render_figure(sim_data,fig_data)
% This function takes data stored externally from sim_save and plots
% it using different settings for different figures. Set up the figure
% settings using the scheme below, and just call this function with the
% figure number and desired data type
%
% Coded by Jeremiah Rhys Wimer, 3/24/2024
close all;

% Set whether the figures are for publication or not

% Print message
fprintf("Generating figure...\n")

% Import data from objects
systems = sim_data.systems;
x_var = sim_data.x_var;
x_range = sim_data.x_range;
num_systems = sim_data.num_systems;
render_for_paper = fig_data.render_for_paper;
save_figure = fig_data.save_figure;
fig_title = fig_data.fig_title;
var_list = fig_data.var_list;
fig_num = fig_data.fig_num;
data_type = fig_data.fig_data_type;
xlabel_vec = fig_data.xlabel_vec;
ylim_vec = fig_data.ylim_vec;
legend_color_vec = fig_data.legend_color_vec;
legend_name_vec = fig_data.legend_name_vec;
final_time = sim_data.final_time;

% Sim loop
results = zeros(num_systems,length(x_range));
for i = 1:num_systems
    % Select current system
    sim_data.current_system = i;
    sys = systems(i);

    % Sweep through all SNR values
    for j = 1:length(x_range)% Select current normalized SNR
        x_val = x_range(j);
        eval("sys." + x_var + " = x_val;");

        % Generate name of data to try to load
        for k = 1:length(var_list)
            if k > 1
                save_name = eval(sprintf("save_name + ""_"" + sys.%s",var_list(k)));
            else
                save_name = eval(sprintf("sys.%s;",var_list(k)));
            end
        end
        load_loc = "Saved Data\\" + save_name + "_sim.mat";

        % Try to load if data exists
        result_file = load(load_loc);
        results_data = result_file.results_data;
        eval(sprintf("results(i,j) = results_data.%s;",data_type));
    end
end

% Make legend command
legend_command = "legend(";
for k = 1:length(legend_name_vec)
    legend_command = legend_command + """" +legend_name_vec(k) + """, ";
end
% legend_command = legend_command + """Location"", ""southwest"");";
legend_command = legend_command + """Location"", """" + fig_data.legend_loc +"""");";

% Change y-axis label to reflect data shown
if data_type == "Thr"
    ylabel_vec = "Throughput (bps/Hz)";
else
    ylabel_vec = data_type;
end

% Plot data
figure(fig_num)
hold on;
for i = 1:num_systems
    try
        % Plot data
        plot(x_range,results(i,:),legend_color_vec(i));
    catch
        % Plot data
        plot(x_range,results(i,:),"-kv");
    end
end
if not(render_for_paper)
    title(fig_title);
else
    line_val = 2;
    mark_val = 10;
    font_val = 16;

    % Adjust plotted line and border line size
    lines = findall(gcf, 'Type', 'Line');
    set(lines, 'LineWidth', line_val);
    ax = gca;
    ax.LineWidth = line_val;

    % Set marker size
    h = findobj(gca, 'Type', 'Line');
    set(h, 'MarkerSize', mark_val);

    % Set font size
    set(ax, 'FontSize', font_val);
end
xlabel(xlabel_vec);
xlim([min(x_range) max(x_range)])
ylabel(ylabel_vec);
ylim(ylim_vec);
if data_type ~= "Thr"
    set(gca, 'YScale', 'log');
end
grid on;
% xticks([3,6,9,12,15,18])
hold off;
eval(legend_command);

if save_figure
    % Save figure with unique name
    if ~exist("Saved Figures", 'dir')
        mkdir("Saved Figures")
    end
    if ~render_for_paper
        saveas(figure(fig_num), sprintf('Saved Figures\\Figure%d_%s.png',fig_num,string(datetime('now', 'format', 'MM.dd.uuuu_HH.mm'))));
    else
        saveas(figure(fig_num), sprintf('Saved Figures\\Figure%d_%s.eps',fig_num,string(datetime('now', 'format', 'MM.dd.uuuu_HH.mm'))), 'epsc');
    end
end

% Completion message
fprintf(sprintf("Figure constructed at %s \n",final_time));