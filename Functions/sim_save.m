function sim_data = sim_save(sim_data)
% This function takes a collection of systems and a SNR range, simulating
% the BER, SER, FER and Throughput of all the systems. The results are
% saved externally to a local folder titled "Saved Data". This
% externally saved data can later be rendered in a figure using the
% render_figure function.
% NOTE: Custom object comms_obj_OTFS is needed for this function to run properly
%
% Coded by Jeremiah Rhys Wimer, 3/24/2024

% Import data from sim_data_obj
systems = sim_data.systems;
num_systems = sim_data.num_systems;
x_range = sim_data.x_range;
max_frames = sim_data.max_frames;
inc_frames = sim_data.inc_frames;
max_trials = sim_data.max_trials;

% Prepare to save data
if ~exist("Saved Data", 'dir')
    mkdir("Saved Data")
end

% Check data in every file for missing data
errors_bit = cell(num_systems,length(x_range));
errors_sym = cell(num_systems,length(x_range));
errors_frame = cell(num_systems,length(x_range));
for i = 1:num_systems
    for j = 1:length(x_range)
        % Select current system
        sys = systems(i,j);
        save_loc = sys.save_loc;

        % Load number of frames and trials
        if isfile(save_loc)
            load_test = load(save_loc);
            results_loaded = load_test.results_data;
            data_bit = results_loaded.bit_errors;
            data_sym = results_loaded.sym_errors;
            data_frame = results_loaded.frame_errors;
        else
            data_bit = NaN;
            data_sym = NaN;
            data_frame = NaN;
        end

        % Adjust size of data as necessary
        [num_frames, num_trials] = size(data_bit);
        maxRows = max(num_frames, max_frames);
        maxCols = max(num_trials, max_trials);
        resized_data = NaN(maxRows, maxCols);

        % Rescale bit array
        resized_data_copy = resized_data;
        resized_data_copy(1:num_frames, 1:num_trials) = data_bit(1:num_frames, 1:num_trials);
        errors_bit{i,j} = resized_data_copy;

        % Rescale sym array
        resized_data_copy = resized_data;
        resized_data_copy(1:num_frames, 1:num_trials) = data_sym(1:num_frames, 1:num_trials);
        errors_sym{i,j} = resized_data_copy;

        % Rescale frame array
        resized_data_copy = resized_data;
        resized_data_copy(1:num_frames, 1:num_trials) = data_frame(1:num_frames, 1:num_trials);
        errors_frame{i,j} = resized_data_copy;
    end
end

% Sim loop
for p = 1:max_trials
    % Set current trial number
    sim_data.current_trial = p;

    for i = 1:num_systems*length(x_range)
        % Print message
        clc;
        fprintf("Checking saved simulation data...\n")

        % Select current system
        sim_data.current_system = i;
        index1 = floor((i-1)/length(x_range))+1;
        index2 = mod(i-1,length(x_range))+1;
        sys = systems(index1,index2);
        data_bit = errors_bit{index1,index2};
        data_sym = errors_sym{index1,index2};
        data_frame = errors_frame{index1,index2};

        % Incrementally render a trial
        sim_data.start_time = datetime('now');
        for j = 1:ceil(max_frames/inc_frames)
            % Set increment range
            frame_range = 1+(j-1)*inc_frames:j*inc_frames;
            if j == ceil(max_frames/inc_frames)
                frame_range = frame_range(frame_range <= max_frames);
            end

            % Check if sim is needed
            sim_flag = false;
            if any(isnan(data_bit(frame_range,p)))
                sim_flag = true;
                sim_data.frame_start = find(isnan(data_bit(frame_range,p)), 1);
                sim_data.frame_range = frame_range;

                sim_data.data_bit = data_bit(frame_range,p);
                sim_data.data_sym = data_sym(frame_range,p);
                sim_data.data_frame = data_frame(frame_range,p);

                end_index = find(isnan(data_bit(:,p)),1) - 1;
                sim_data.carry_bit = sum(data_bit(1:end_index,p),"all");
                sim_data.carry_sym = sum(data_sym(1:end_index,p),"all");
                sim_data.carry_frame = sum(data_frame(1:end_index,p),"all");
            end

            % Run simulation
            if sim_flag
                % Simulation function
                [sys,sim_data] = OTFS_pulse_sim(sys,sim_data);
                systems(index1,index2) = sys;

                % Update local data
                data_bit(frame_range,p) = sim_data.data_bit;
                data_sym(frame_range,p) = sim_data.data_sym;
                data_frame(frame_range,p) = sim_data.data_frame;
                errors_bit{index1,index2} = data_bit;
                errors_sym{index1,index2} = data_sym;
                errors_frame{index1,index2} = data_frame;
            end

            if sim_flag || p==max_trials
                % Save data
                results_data = obj_results;
                results_data.max_frames = max_frames;
                results_data.max_trials = max_trials;
                results_data.M_ary = sys.M_ary;
                results_data.syms_per_frame = sys.syms_per_f;
                results_data.bit_errors = data_bit;
                results_data.sym_errors = data_sym;
                results_data.frame_errors = data_frame;
                save(sys.save_loc,"results_data")
            end
        end
    end
end

% Announce completion time
sim_data.final_time = datetime('now');