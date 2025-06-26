classdef obj_sim
    %% SIM_OBJ
    %% How to use
    % This object stores the simulation information outside the simulated
    % systems data, so it can update needed parameters like needed runtime
    % and simulation display information.
    %
    % Written by JRW, 6/21/2024

    %% Properties ---------------------------------------------------------
    properties
        % Display parameters
        freq_display;
        num_systems = 1;
        current_system = 0;
        current_trial = 0;
        frame_start;
        start_time;
        
        % Moving data in and out of sim
        data_bit;
        data_sym;
        data_frame;
        frame_range;
        carry_bit;
        carry_sym;
        carry_frame;

        % Simulation limits for each system
        test_select;
        max_frames = 250;
        inc_frames = 50;
        max_trials = 10;
        tol = 0;
        signif_level = 0.05;

        % Simulation parameters
        need_sim = true;
        systems;
        var_list;
        var1 = "N_tsyms";
        var2 = "NA";
        x_var = "EbN0_db";
        x_range = 3:3:18;
        var1_range = [16,64];
        var2_range = ["rect","sinc","rrc"];
        var_defaults = {18,4,"MPSK",15e3,4e9,500,16,64,8,"rect",1};
        rrc_vals = 1;
        init_time;
        final_time;

        % Reload settings
        reload = true;
        % save_loc;
    end

    properties (Dependent)
        total_sims
    end

    methods
        %% DEPENDENT VARIABLES --------------------------------------------

        function result = get.total_sims(obj)
            result = length(obj.x_range) * obj.num_systems;
        end

        function obj = sim_obj()
            obj.init_time = datetime('now');
        end

        %% INTERNAL FUNCTIONS ---------------------------------------------

        
    end
end