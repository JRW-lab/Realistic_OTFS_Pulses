%% Initialization
clear; clc; close all;
addpath(fullfile(pwd, 'Functions'));

% Initialization
sim_data = obj_sim;
fig_data = obj_fig;
profile_num = 0;

% System setting defaults - list all possible variables
var_list = ["EbN0_db","M_ary","select_mod","sbcar_spacing","Fc","v_vel","N_tsyms","M_sbcars","q","filter","rolloff"];
fig_data.var_list = var_list;

%% Simulation and figure settings
% Simulation settings
freq_of_display = 1;        % Hz console updates
max_frames = 100;           % Max frames per trial
max_trials = 100;            % Max trials per system
increment_frames = 25;      % Saves at mod(frames,increment_frames) and the end
reload = true;              % Reload previously-saved data

% Figure settings
fig_num = 1;
make_all_figs = false;
show_figure = true;
save_figure = true;
render_for_paper = true;

% Figure 1 - BER vs all filters, Q=8
profile_num = profile_num + 1;
new_profile = obj_profile;
sub_sim_data = sim_data;
sub_fig_data = fig_data;
sub_sim_data.var1 = "N_tsyms";
sub_sim_data.var2 = "filter";
sub_sim_data.x_var = "EbN0_db";
sub_sim_data.x_range = 3:3:18;
sub_sim_data.var1_range = [16,64];
sub_sim_data.var2_range = ["rect","sinc","RRC","RRC"];
sub_sim_data.var_defaults = {18,4,"MPSK",15e3,4e9,500,16,64,8,"rect",1};
sub_sim_data.rrc_vals = [1,0.2];
sub_fig_data.fig_num = profile_num;
sub_fig_data.fig_data_type = "BER";
sub_fig_data.xlabel_vec = "E_b/N_0 (dB)";
sub_fig_data.ylim_vec = [1e-7 1e-1];
sub_fig_data.fig_title = "Performance of typical pulses under Q=8";
sub_fig_data.var1_spec = ["r","b"];
sub_fig_data.var2_spec = ["-+","--o","-.*",":v"];
sub_fig_data.gen_spec = "";
new_profile.sub_sim_data = sub_sim_data;
new_profile.sub_fig_data = sub_fig_data;
eval("profile_"+profile_num+" = new_profile;");

% Figure 2 - BER vs Q, sinc and RRC
profile_num = profile_num + 1;
new_profile = obj_profile;
sub_sim_data = sim_data;
sub_fig_data = fig_data;
sub_sim_data.var1 = "rolloff";
sub_sim_data.var2 = "NA";
sub_sim_data.x_var = "q";
sub_sim_data.x_range = 2:2:14;
sub_sim_data.var1_range = [0,0.2,0.5,0.8,1];
sub_sim_data.var2_range = [];
sub_sim_data.var_defaults = {15,4,"MPSK",15e3,4e9,500,16,64,8,"rrc",1};
sub_sim_data.rrc_vals = sub_sim_data.var1_range;
sub_fig_data.fig_num = profile_num;
sub_fig_data.fig_data_type = "BER";
sub_fig_data.xlabel_vec = "Q";
sub_fig_data.ylim_vec = [1e-4 2e-2];
sub_fig_data.fig_title = "Improvements with RRC pulses as Q increases";
sub_fig_data.var1_spec = ["-r+","--ro","-b*","--bv","-ksquare"];
sub_fig_data.var2_spec = "";
sub_fig_data.gen_spec = "";
sub_fig_data.legend_loc = "northeast";
new_profile.sub_sim_data = sub_sim_data;
new_profile.sub_fig_data = sub_fig_data;
eval("profile_"+profile_num+" = new_profile;");

% Figure 3 - BER vs rolloff, RRC
profile_num = profile_num + 1;
new_profile = obj_profile;
sub_sim_data = sim_data;
sub_fig_data = fig_data;
sub_sim_data.var1 = "N_tsyms";
sub_sim_data.var2 = "NA";
sub_sim_data.x_var = "rolloff";
sub_sim_data.x_range = 0:.1:1;
sub_sim_data.var1_range = [16,32,64];
sub_sim_data.var2_range = [];
sub_sim_data.var_defaults = {15,4,"MPSK",15e3,4e9,500,16,64,8,"rrc",1};
sub_fig_data.fig_num = profile_num;
sub_fig_data.fig_data_type = "BER";
sub_fig_data.xlabel_vec = "rolloff factor \alpha";
sub_fig_data.ylim_vec = [3e-5 1e-3];
sub_fig_data.fig_title = "Performance of RRC under various \alpha values";
sub_fig_data.var1_spec = ["-r+","--ko","-.b*"];
sub_fig_data.var2_spec = "";
sub_fig_data.gen_spec = "";
sub_fig_data.legend_loc = "northwest";
new_profile.sub_sim_data = sub_sim_data;
new_profile.sub_fig_data = sub_fig_data;
eval("profile_"+profile_num+" = new_profile;");

% Figure 4 - Resiliency to velocity
profile_num = profile_num + 1;
new_profile = obj_profile;
sub_sim_data = sim_data;
sub_fig_data = fig_data;
sub_sim_data.var1 = "filter";
sub_sim_data.var2 = "NA";
sub_sim_data.x_var = "v_vel";
sub_sim_data.x_range = 50:150:800;
sub_sim_data.var1_range = ["rect","sinc","RRC","RRC"];
sub_sim_data.var2_range = [];
sub_sim_data.var_defaults = {15,4,"MPSK",15e3,4e9,500,16,64,8,"rect",0.3};
sub_sim_data.rrc_vals = [1,0.2];
sub_fig_data.fig_num = profile_num;
sub_fig_data.fig_data_type = "BER";
sub_fig_data.xlabel_vec = "vehicle speed (km/hr)";
sub_fig_data.ylim_vec = [1e-5 1e-2];
sub_fig_data.fig_title = "Performance under various vehicle speeds";
sub_fig_data.var1_spec = ["r-+","b-o","k-*","k-v"];
sub_fig_data.var2_spec = "";
sub_fig_data.gen_spec = "";
new_profile.sub_sim_data = sub_sim_data;
new_profile.sub_fig_data = sub_fig_data;
eval("profile_"+profile_num+" = new_profile;");

%% Backend
% Set figures that are simulated and rendered
fprintf("Selecting figures to be rendered...\n")
if not(make_all_figs)
    sweep_range = fig_num;
else
    sweep_range = 1:profile_num;
end

% Set up data for all profiles
fprintf("Creating profiles from given user inputs...\n")
all_profiles = repmat(obj_profile,profile_num,1);
for i = 1:profile_num
    eval("all_profiles(i) = profile_" + i + ";");
end
for p = sweep_range
    % Set current profile
    profile_sel = all_profiles(p);

    % Set sim data
    fprintf(sprintf("   Setting sim data for Profile %d...\n",p))
    sim_data = profile_sel.sub_sim_data;
    sim_data.max_trials = max_trials;
    sim_data.max_frames = max_frames;
    sim_data.inc_frames = increment_frames;
    sim_data.freq_display = freq_of_display;
    sim_data.var_list = var_list;
    sim_data.reload = reload;

    % Set fig data
    fprintf(sprintf("   Setting fig data for Profile %d...\n",p))
    fig_data = profile_sel.sub_fig_data;
    fig_data.render_for_paper = render_for_paper;
    fig_data.save_figure = save_figure;

    % Extract sim data to make systems
    fprintf(sprintf("   Creating parameters for Profile %d...\n",p))
    var1 = sim_data.var1;
    var2 = sim_data.var2;
    x_var = sim_data.x_var;
    x_range = sim_data.x_range;
    var1_range = sim_data.var1_range;
    var2_range = sim_data.var2_range;
    var_defaults = sim_data.var_defaults;

    % Make needed variables
    if var2 == "NA"
        num_vars = 1;
        var_names = var1;
        variants = table2cell(combinations(var1_range));
    else
        num_vars = 2;
        var_names = [var1,var2];
        variants = table2cell(combinations(var1_range,var2_range));
    end

    % Add a detail about selected alpha for RRC pulses
    rrc_count = 0;
    for k = 1:2
        try
            vars_temp = [variants{:,k}];
            if ismember("rrc",vars_temp) || ismember("RRC",vars_temp)
                for i = 1:length(vars_temp)
                    if vars_temp(i) == "rrc" || vars_temp(i) == "RRC"
                        rrc_index = mod(rrc_count,length(sim_data.rrc_vals)) + 1;
                        vars_temp(i) = vars_temp(i) + ", \alpha=" + sim_data.rrc_vals(rrc_index);
                        rrc_count = rrc_count + 1;
                    end
                end
            end

            if k == 1
                if size(variants,2) == 1
                    variants_legend = table2cell(combinations(vars_temp));
                else
                    variants_legend = table2cell(combinations(vars_temp(1:(length(vars_temp)/length(var2_range))),var2_range));
                end
            else
                variants_legend = table2cell(combinations(var1_range,vars_temp(1:(length(vars_temp)/length(var1_range)))));
            end

            break;
        catch
            variants_legend = variants;
        end
    end

    % Make legend vector
    fprintf(sprintf("   Creating legend vector for Profile %d...\n",p))
    legend_vec_name = repmat("",1,size(variants_legend,1));
    legend_mat = cell(size(variants_legend,1),num_vars);
    for j = 1:num_vars
        var_sel = eval("var" + j);
        switch var_sel
            case "EbN0_db"
                prefix = "";
                suffix = "dB";
            case "M_ary"
                prefix = "";
                suffix = "-ary";
            case "select_mod"
                prefix = "";
                suffix = "";
            case "sbcar_spacing"
                prefix = "df=";
                suffix = "Hz";
            case "Fc"
                prefix = "f_c=";
                suffix = "Hz";
            case "v_vel"
                prefix = "";
                suffix = "km/hr";
            case "N_tsyms"
                prefix = "N=";
                suffix = "";
            case "M_sbcars"
                prefix = "M=";
                suffix = "";
            case "q"
                prefix = "Q=";
                suffix = "";
            case "filter"
                prefix = "";
                suffix = "";
            case "rolloff"
                prefix = "\alpha=";
                suffix = "";
        end
        for i = 1:size(variants_legend,1)
            legend_mat{i,j} = prefix + variants_legend{i,j} + suffix;
        end
    end
    for i = 1:size(variants_legend,1)
        for j = 1:num_vars
            if j > 1
                legend_vec_name(i) = legend_vec_name(i) + ",";
            end
            legend_vec_name(i) = legend_vec_name(i) + legend_mat{i,j};
        end
    end
    fig_data.legend_name_vec = legend_vec_name;

    % Define all systems
    fprintf(sprintf("   Creating system objects for Profile %d...\n",p))
    num_variants = size(variants,1);
    systems = repmat(obj_comms_OTFS,num_variants,1);
    for i = 1:length(systems)
        for k = 1:length(var_list)
            if ismember(var_list(k),var_names)
                variant_index = find(strcmp(var_list(k),var_names));
                variant_val = variants(i,variant_index);
                eval("systems(i)." + var_list(k) + " = variant_val{1};");
            else
                % variant_val = var_defaults{k};
                eval("systems(i)." + var_list(k) + " = var_defaults{k};");
            end
        end
    end

    % Add redundancies
    fprintf(sprintf("   Creating system objects for Profile %d...\n",p))
    rrc_count = 0;
    for i = 1:length(systems)
        if systems(i).filter == "rect"
            systems(i).q = 1;
            systems(i).rolloff = 1;
        end
        if systems(i).filter == "sinc"
            systems(i).rolloff = 1;
        end
        if systems(i).filter == "rrc" || systems(i).filter == "RRC"
            if systems(i).rolloff == 0
                systems(i).filter = "sinc";
                systems(i).rolloff = 1;
                rrc_count = rrc_count + 1;
            else
                rrc_index = mod(rrc_count,length(sim_data.rrc_vals)) + 1;
                systems(i).rolloff = sim_data.rrc_vals(rrc_index);
                rrc_count = rrc_count + 1;
            end
        end
    end

    % Make new variants of systems based on x axis and assign filename
    fprintf(sprintf("   Creating parameter sweep for Profile %d...\n",p))
    systems = repmat(systems,1,length(x_range));
    for i = 1:size(systems,1)
        for j = 1:size(systems,2)
            % Select current x-val to modify system
            x_val = x_range(j);
            eval("systems(i,j)." + x_var + " = " + x_val + ";");

            % Generate name of data to try to load
            for k = 1:length(var_list)
                if k > 1
                    save_name = eval(sprintf("save_name + ""_"" + systems(i,j).%s",var_list(k)));
                else
                    save_name = eval(sprintf("systems(i,j).%s;",var_list(k)));
                end
            end
            save_loc = "Saved Data\\" + save_name + "_sim.mat";
            systems(i,j).save_loc = save_loc;

        end
    end

    % Assign systems back to sim data
    fprintf(sprintf("   Assigning system data back to Profile %d...\n",p))
    num_systems = size(systems,1);
    sim_data.num_systems = num_systems;
    sim_data.systems = systems;
    
    % Assign metadata back to profiles
    fprintf(sprintf("   Assigning metadata back to Profile %d...\n",p))
    profile_sel.sub_sim_data = sim_data;
    profile_sel.sub_fig_data = fig_data;
    all_profiles(p) = profile_sel;
end

% Sweep through all systems
for p = sweep_range
    % Extract profile and needed data
    profile_sel = all_profiles(p);
    sim_data = profile_sel.sub_sim_data;
    fig_data = profile_sel.sub_fig_data;

    % Run simulations as needed
    sim_data = sim_save(sim_data);

    % Render figure if end is shown
    if show_figure
        % Plot data
        render_figure(sim_data,fig_data);
    end
end
