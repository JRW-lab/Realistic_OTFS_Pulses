function [sys,sim_data] = OTFS_pulse_sim(sys,sim_data)
% This is the actual simulation file, configured for the shuffled
% doppler-Domain OTFS system model with truncated pulse shapes
% (rect,sinc,rrc)
%
% Written by JRW, 6/21/2024

% % Debug inputs
% clear; clc;
% obj = comms_obj_OTFS;
% obj.Eb_N0_db = 5;
% obj.M_ary = 4;
% obj.N_tsyms = 16;
% obj.M_sbcars= 64;
% obj.filter = "rrc";
% obj.q = 2;
% max_errors = 1000;
% max_frames = 5000;
% max_time = Inf;
% settings_number = 1;

% Equalizer settings
N_iters = 8;
ambig_res = 1001;
% rng('default');

% Inputs from sim_data object
num_systems = sim_data.num_systems * length(sim_data.x_range);
current_system = sim_data.current_system;
current_trial = sim_data.current_trial;
freq_display = sim_data.freq_display;
max_trials = sim_data.max_trials;
max_frames = sim_data.max_frames;
bit_errors = sim_data.data_bit;
sym_errors = sim_data.data_sym;
frame_errors = sim_data.data_frame;
frame_range = sim_data.frame_range;
frame_start = frame_range(1);
carry_bit = sim_data.carry_bit;
carry_sym = sim_data.carry_sym;
carry_frame = sim_data.carry_frame;
start_time = sim_data.start_time;

% Inputs from system object
EbN0 = sys.EbN0_db;
Es = sys.Es;
S = sys.S;
M_ary = sys.M_ary;
N0 = sys.N0;
% N0 = 0; % Uncomment to check noiseless scenario
M = sys.M_sbcars;
N = sys.N_tsyms;
Fc = sys.Fc;
v = sys.v_vel;
if M_ary == 2
    bit_order = [0;1];
elseif M_ary == 4
    bit_order = [0,0;0,1;1,0;1,1];
end
shape = sys.filter;
alpha = sys.rolloff;
q = sys.q;
T = sys.T;
Ts = sys.Ts;
Lp = sys.Lp;
Ln = sys.Ln;
if shape == "rect"
    q = 1;
end

load_test = sys.ambig_vals;
if isempty(load_test)
    % Set discrete ambiguity lookup table if not done already
    clc;
    fprintf(sprintf("Checking ambiguity table and noise covariance matrix for system %d of %d...\n",current_system,num_systems))
    if shape ~= "rect"
        if ~exist("Saved Data\\Discrete Ambiguity Tables", 'dir')
            mkdir("Saved Data\\Discrete Ambiguity Tables")
        end
        filename = sprintf("Saved Data\\Discrete Ambiguity Tables\\ambig_discrete_q%d_%s_alpha%.1f.mat",q,shape,alpha);
        if isfile(filename)
            loaded_file = load(filename);
            ambig_t_range = loaded_file.ambig_t_range;
            ambig_f_range = loaded_file.ambig_f_range;
            ambig_vals = loaded_file.ambig_vals;
        else
            ambig_t_lim = Ts*(q + ceil((2510*10^(-9))/Ts));
            ambig_f_lim = ((v * (1000/3600))*Fc) / (physconst('LightSpeed'));
            ambig_t_range = linspace(-ambig_t_lim,ambig_t_lim,ambig_res);
            ambig_f_range = linspace(-ambig_f_lim,ambig_f_lim,ambig_res);
            ambig_vals = zeros(ambig_res);
            update_vals = floor(length(ambig_t_range)*(linspace(.01,.99,34)));
            clc;
            fprintf("Generating discrete ambiguity table (one time process)...\nProgress:                        |\n")
            for k = 1:length(ambig_t_range)
                if ismember(k,update_vals)
                    fprintf("x")
                end

                for l = 1:length(ambig_f_range)
                    ambig_vals(k,l) = ambig_direct(ambig_t_range(k),ambig_f_range(l),Ts,shape,alpha,q,ambig_res);
                end
            end
            fprintf("\n")

            % Save to file
            save(filename,"ambig_t_range","ambig_f_range","ambig_vals");
        end
    else
        ambig_t_range = [];
        ambig_f_range = [];
        ambig_vals = [];
    end
    sys.ambig_t_range = ambig_t_range;
    sys.ambig_f_range = ambig_f_range;
    sys.ambig_vals = ambig_vals;
end

load_test = sys.R;
if isempty(load_test)
    % Set up noise covariance matrix if not done already
    if shape ~= "rect"
        if ~exist("Saved Data\\Noise Covariance Matrices", 'dir')
            mkdir("Saved Data\\Noise Covariance Matrices")
        end
        filename = sprintf("Saved Data\\Noise Covariance Matrices\\Rzddt_N%d_M%d_q%d_%s_alpha%.1f.mat",N,M,q,shape,alpha);
        if isfile(filename)
            loaded_file = load(filename);
            R_half = loaded_file.R_half;
            R = loaded_file.R;
        else
            clc;
            fprintf("Generating noise covariance matrix (one time process)...\n")
            [R_half,R] = gen_Rzddt(T,N,M,q,ambig_res,shape,alpha);
            fprintf("\n")

            % Save to file
            save(filename,"R_half","R");
        end
    else
        R_half = eye(N*M);
        R = R_half;
    end
    sys.R_half = R_half;
    sys.R = R;
end

% Needed variables and matrices setup
syms_per_f = M*N;
Gamma_MN = gen_Gamma_MN(M,N);
F_N = gen_DFT(N);

% Bring in previously loaded data
ambig_t_range = sys.ambig_t_range;
ambig_f_range = sys.ambig_f_range;
ambig_vals = sys.ambig_vals;
R = sys.R;
R_half = sys.R_half;
if shape ~= "rect"
    res_chn_tau = (ambig_t_range(2)-ambig_t_range(1));
    res_chn_v = ambig_f_range(2)-ambig_f_range(1);
end

% Start sim loop
timer = tic;
display_flag = false;
frame_display = Inf;

% Simulation loop
for frames = frame_range

    % Generate data
    [TX_1,TX_2,x_DD] = gen_data(bit_order,S,syms_per_f);
    TX_bit = Gamma_MN' * TX_1;
    TX_sym = Gamma_MN' * TX_2;
    x_tilde = Gamma_MN' * x_DD;

    % Generate channel
    [chn_g,chn_tau,chn_v] = channel_generation(Fc,v);
    if shape == "rect" % rectangular ambiguity is closed form
        % Create H Matrix
        H = gen_H(T,N,M,Lp,Ln,chn_g,chn_tau,chn_v,shape,alpha);
    else
        % Normalize tau and v to cohere with discrete ambig values
        chn_tau = round(chn_tau/res_chn_tau)*res_chn_tau;
        chn_v = round(chn_v/res_chn_v)*res_chn_v;

        % Find direct tap indices and tap values
        l = (Ln:Lp).';
        tap_t_range = (l*Ts - chn_tau) .* ones(Lp-Ln+1,length(chn_g),N*M);
        tap_f_range = (ones(Lp-Ln+1,1) .* chn_v) .* ones(Lp-Ln+1,length(chn_g),N*M);
        tap_t_range = round(tap_t_range ./ res_chn_tau) + ceil(length(ambig_t_range)/2);
        tap_f_range = round(tap_f_range ./ res_chn_v) + ceil(length(ambig_f_range)/2);
        tap_t_range(tap_t_range < 1) = 1;
        tap_t_range(tap_t_range > 1001) = 1001;
        tap_f_range(tap_f_range < 1) = 1;
        tap_f_range(tap_f_range > 1001) = 1001;

        % Create H Matrix
        H = gen_H_direct(T,N,M,Lp,Ln,chn_g,chn_v,ambig_vals,tap_t_range,tap_f_range);
    end
    H_tilde = gen_H_tilde(H,M,N,Lp,Ln,F_N);

    % Generate noise
    z_tilde = sqrt(N0/2) * R_half * (randn(syms_per_f,1) + 1j*randn(syms_per_f,1));

    % Create receive vector
    y_tilde = H_tilde * x_tilde + z_tilde;

    % Iterative Detector - JRW
    % x_hat = OTFS_pulse_equalizer(y_tilde,H_tilde,N,M,Lp,Ln,Es,N0,S,N_iters);
    x_hat = OTFS_pulse_equalizer_v2(y_tilde,H_tilde,N,M,Lp,Ln,Es,N0,S,N_iters,R);

    % Hard detection for final x_hat
    dist = abs(x_hat.' - S).^2;
    [~,min_index] = min(dist);
    RX_sym = min_index.' - 1;

    % Convert final RX_sym to RX_bit
    RX_bit = bit_order(RX_sym+1,:);

    % Error calculation
    bit_error_vec = TX_bit ~= RX_bit;
    sym_error_vec = TX_sym ~= RX_sym;
    bit_errors(frames-frame_start+1) = sum(bit_error_vec,"all");
    sym_errors(frames-frame_start+1) = sum(sym_error_vec,"all");
    if sum(bit_error_vec,"all") > 0
        frame_errors(frames-frame_start+1) = 1;
    else
        frame_errors(frames-frame_start+1) = 0;
    end

    % Record sim time to set display update frequency
    run_time = toc(timer);
    sim_time = run_time/(frames-frame_start+1);
    if run_time > freq_display && not(display_flag)
        frame_display = frames-frame_start+1;
        display_flag = true;
    end

    % Update progress status
    if (mod((frames-frame_start+1),frame_display) == 0 && display_flag) || frames == frame_range(end)
        % Set rates when time for it
        bit_errors_so_far = sum(bit_errors(1:(frames-frame_start+1)),"all") + carry_bit;
        sym_errors_so_far = sum(sym_errors(1:(frames-frame_start+1)),"all") + carry_sym;
        frame_errors_so_far = sum(frame_errors(1:(frames-frame_start+1)),"all") + carry_frame;
        BER = bit_errors_so_far / (frames * syms_per_f * log2(M_ary));
        SER = sym_errors_so_far / (frames * syms_per_f);
        FER = frame_errors_so_far / frames;

        % Get needed variables
        projected_length = sim_time * max_frames;
        projected_end_time = start_time + seconds(projected_length);

        % Print message
        clc;
        fprintf(sprintf("Simulation initialized: %s \n", char(start_time)))
        fprintf(sprintf("Projected end time:     %s \n", char(projected_end_time)))
        fprintf(sprintf("Trial %d of %d > System %d of %d > Frame %d of %d \n", current_trial, max_trials, current_system, num_systems, frames, max_frames))
        fprintf(sprintf("Simulation parameters: %.2fdB, %dPSK, df=%.0f, M=%d, N=%d, v=%d, %s, Q=%d, alpha=%.1f \n", EbN0,M_ary,1/T,M,N,v,shape,q,alpha))
        fprintf(sprintf("        BER: %.3d | %d bit errors\n", BER, bit_errors_so_far))
        fprintf(sprintf("        SER: %.3d | %d symbol errors\n", SER, sym_errors_so_far))
        fprintf(sprintf("        FER: %.3d | %d frame errors\n", FER, frame_errors_so_far))
    end
end

% Overwrite sim data with new frame info
sim_data.data_bit = bit_errors;
sim_data.data_sym = sym_errors;
sim_data.data_frame = frame_errors;
