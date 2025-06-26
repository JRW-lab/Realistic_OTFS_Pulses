%% MAIN_CONTINUOUS_SAMPLE
% This MATLAB file simulates doppler-delay domain information, produces a
% continuous signal passing through a rect/sinc/RRC filter, passes through
% a doubly-selective (t and f) channel, passes through another filter then
% samples.
% The resulting sampled RX vector is then directly compared with the
% simpler matrix-based method.
%
% Written by JRW, 5/20/2024

clc; close all;

% System settings
show_figures = false;
res_discretization = 10;
q = 2;
alpha = 1;
num_iters = 10;
Eb_N0_db = Inf;
M = 64;
N = 16;

% System initialization
Eb = 1;
N0 = (Eb / (10^(Eb_N0_db / 10)));
bit_order = [0,0;0,1;1,0;1,1];
S = Eb * [(sqrt(2)/2) + (1j*sqrt(2)/2);...
     (sqrt(2)/2) - (1j*sqrt(2)/2);...
     -(sqrt(2)/2) + (1j*sqrt(2)/2);...
     -(sqrt(2)/2) - (1j*sqrt(2)/2)];
T = 1/15000;
Ts = T/M;
dt = Ts / res_discretization;
syms_per_f = M*N;
Fc = 4*10^9;
v = 500;
pulse_names = ["rect","sinc","rrc"].';

% Generate channel info for Lp-Ln taps and normalize noise
Lp =  q + floor((2510*10^(-9)) / Ts);
Ln = -q;
N0 = N0 * (N*M+Lp-Ln) / (N*M);

% Set up block transforms
Gamma_MN = gen_Gamma_MN(M,N);
F_M = gen_DFT(M);
F_N = gen_DFT(N);

% Set up time ranges
t_range_sig = -Ts*Lp:dt:N*T-Ts*Ln-dt;
t_range_filter = 0:dt:q*Ts;
t_range_tx = -Ts*(Lp):dt:N*T+Ts*(-Ln+q)-dt;
t_range_rx = -Ts*(Lp+q):dt:N*T+Ts*(-Ln+q)-dt;

%% Filters
fprintf("Generating rect/sinc/rrc filters...\n")
% Filter setup
filters = zeros(3,length(t_range_filter));
filters(1,:) = gen_rect_filter(t_range_filter,Ts,q/2);
filters(2,:) = sinc_trunc((t_range_filter - q*Ts/2),Ts,q/2);
filters(3,:) = RRCt_trunc((t_range_filter - q*Ts/2),alpha,Ts,q/2);

% Normalize E of filters
for k = 1:3
    filters(k,:) = filters(k,:) ./ sqrt(sum(filters(k,:).^2));
end

%% Discretize Channel Info
fprintf("Generating ambiguity value tables for use later...\n")
% Set up ambiguity lookup table
max_Delay = Ts*(q + ceil((2510*10^(-9))/Ts));
max_Doppler = ((v * (1000/3600))*Fc) / (physconst('LightSpeed'));
ambig_t_range = (-max_Delay:dt:max_Delay);
ambig_f_range = linspace(-max_Doppler,max_Doppler,length(ambig_t_range));
ambig_vals = zeros(length(ambig_t_range),length(ambig_f_range),3);
for k = 1:length(ambig_t_range)
    for l = 1:length(ambig_f_range)
        for m = 1:3
            ambig_vals(k,l,m) = ambig_direct(ambig_t_range(k),ambig_f_range(l),Ts,pulse_names(m),alpha,q,res_discretization);
        end
    end
end

% Generate channel info to use later
[chn_g,chn_tau,chn_v] = channel_generation(Fc,v);

% Define resolution of tau and v
res_chn_tau = (ambig_t_range(2)-ambig_t_range(1));
res_chn_v = ambig_f_range(2)-ambig_f_range(1);

% Normalize tau and v to cohere with discrete ambig values
chn_tau = round(chn_tau/res_chn_tau)*res_chn_tau;
chn_v = round(chn_v/res_chn_v)*res_chn_v;

% Find direct tap indices and values
l = (Ln:Lp).';
tap_t_range = (l*Ts - chn_tau) .* ones(Lp-Ln+1,length(chn_g),N*M);
tap_f_range = (ones(Lp-Ln+1,1) .* chn_v) .* ones(Lp-Ln+1,length(chn_g),N*M);
tap_t_indices = round(tap_t_range ./ res_chn_tau) + ceil(length(ambig_t_range)/2);
tap_f_indices = round(tap_f_range ./ res_chn_v) + ceil(length(ambig_f_range)/2);

%% Data
% Generate data
fprintf("Generating a continuous time-domain signal...\n")
[~,TX_2,x_dd_vec] = gen_data(bit_order,S,syms_per_f);
TX_sym = Gamma_MN' * TX_2;
x_tilde = Gamma_MN' * x_dd_vec;

% Get relevent vectors and matrices from data
x_dd_block = zeros(M,N);
for i = 0:N*M-1
    dim1 = mod(i,M)+1;
    dim2 = floor(i/M)+1;
    x_dd_block(dim1,dim2) = x_dd_vec(i+1);
end
x_tf_block = F_M * x_dd_block * F_N';
sig_block = F_M' * x_tf_block;
sig_vec = zeros(N*M,1);
for i = 0:N*M-1
    dim1 = mod(i,M)+1;
    dim2 = floor(i/M)+1;
    sig_vec(i+1) = sig_block(dim1,dim2);
end

% Add prefix and postfix for convolutions
sig_prefix = sig_vec(N*M-Lp+1:N*M);
sig_postfix = sig_vec(1:-Ln);
sig_buffer = [sig_prefix; sig_vec; sig_postfix];

% Make transmitted signal a train of deltas
sig_t = zeros(numel(sig_buffer)*res_discretization,1);
sig_t(1:res_discretization:end) = sig_buffer;

% View original signal
if show_figures
    figure(1)
    title("TX Signal Before Pulse Shaping")
    hold on
    plot(t_range_sig/Ts,real(sig_t))
    plot(t_range_sig/Ts,imag(sig_t))
    hold off
    grid on;
    legend("Real Part","Imag Part")
end

%% TX
fprintf("Passing data through first filter...\n");
% Convolute time domain TX signal with first filter
s_t = zeros(length(t_range_tx),3);
for k = 1:3
    s_t(:,k) = conv(filters(k,:),sig_t);
end

% View TX signals
if show_figures
    figure(2)
    for i = 1:3
        subplot(3,1,i)
        title_string = sprintf("TX Signal (%s Pulse)",pulse_names(i));
        title(title_string);
        hold on;
        plot(t_range_tx/Ts,real(s_t(:,i)))
        plot(t_range_tx/Ts,imag(s_t(:,i)))
        grid on;
        legend("Real Part","Imag Part")
        hold off;
    end
end

%% Channel
% Pass filtered signal through channel
r_t = zeros(length(t_range_tx),3);
update_vals = floor(length(t_range_tx)*(linspace(.01,.99,34)));
fprintf("Passing differently filtered signals through channel...\nProgress:                        |\n")
for t_index = 1:length(t_range_tx)-1
    if ismember(t_index,update_vals)
        fprintf("x")
    end

    % Take specific t instance
    t = t_range_tx(t_index);

    % Run through all sums of the channel
    rx_sums = zeros(3,1);
    for i = 1:length(chn_v)
        s_time = t - chn_tau(i);
        [~,s_index] = min(abs(t_range_tx - s_time).^2);
        for k = 1:3
            sum_val = chn_g(i) * s_t(s_index,k) * exp(1j*2*pi*chn_v(i)*(t-chn_tau(i)));
            rx_sums(k) = rx_sums(k) + sum_val;
        end
    end
    % Set channel coefficient for each filter
    for k = 1:3
        r_t(t_index,k) = rx_sums(k);
    end
end
fprintf("\nComplete!\n");

% Add noise to RX signal
z_t = sqrt(N0/2) * (rand(length(r_t(:,1)),1) + 1j*rand(length(r_t(:,1)),1));
r_t = r_t + z_t;

%% RX
% Convolute time-domain RX signal and noise with second filter
fprintf("Filtering and sampling received signal...\n")
eta_t = zeros(length(t_range_rx),3);
w_t = eta_t;
for k = 1:3
    eta_t(:,k) = conv(flip(filters(k,:)),r_t(:,k));
    w_t(:,k)  = conv(flip(filters(k,:)),z_t);
end

% View RX signals
if show_figures
    figure(3)
    for i = 1:3
        subplot(3,1,i)
        title_string = sprintf("RX Signal (%s Pulse)",pulse_names(i));
        title(title_string);
        hold on;
        plot(t_range_rx/Ts,real(eta_t(:,i)))
        plot(t_range_rx/Ts,imag(eta_t(:,i)))
        grid on;
        legend("Real Part","Imag Part")
        hold off;
    end
end

% Sample RX signal after filtering for discrete representation
eta_vec_ext = zeros(N*M+Lp-Ln+2*q,3);
w_vec_ext = eta_vec_ext;
for k = 1:3
    eta_vec_ext(:,k) = eta_t(1:res_discretization:end,k);
    w_vec_ext(:,k) = w_t(1:res_discretization:end,k);
end

% Truncate sampled vector to just one frame
eta_cont = zeros(N*M,3);
w_vec = zeros(N*M,3);
for k = 1:3
    eta_cont(:,k) = eta_vec_ext(1+Lp+q:N*M+Lp+q,k);
    w_vec(:,k) = w_vec_ext(1+Lp+q:N*M+Lp+q,k);
end

% Get relevent vectors and matrices from noise
w_block = zeros(M,N,3);
for i = 0:N*M-1
    dim1 = mod(i,M)+1;
    dim2 = floor(i/M)+1;
    for k = 1:3
        w_block(dim1,dim2,k) = w_vec(i+1,k);
    end
end
z_tilde = zeros(N*M,3);
for k = 1:3
    z_tf_block = F_M * w_block(:,:,k);
    z_dd_block = F_M' * z_tf_block * F_N;
    z_dd_vec = zeros(N*M,1);
    for i = 0:N*M-1
        dim1 = mod(i,M)+1;
        dim2 = floor(i/M)+1;
        z_dd_vec(i+1) = z_dd_block(dim1,dim2);
    end
    z_tilde(:,k) = Gamma_MN' * z_dd_vec;
end

% Create block versions of eta for IDFT and SFFT operations
eta_block = zeros(M,N,3);
for i = 0:N*M-1
    dim1 = mod(i,M)+1;
    dim2 = floor(i/M)+1;
    for k = 1:3
        eta_block(dim1,dim2,k) = eta_cont(i+1,k);
    end
end

% Convert to TF domain then DD
y_tf_block = zeros(M,N,3);
y_dd_block = zeros(M,N,3);
for k = 1:3
    y_tf_block(:,:,k) = F_M * eta_block(:,:,k);
    y_dd_block(:,:,k) = F_M' * y_tf_block(:,:,k) * F_N;
end

% Convert back to vectors
y_dd_vec = zeros(N*M,3);
for i = 0:N*M-1
    dim1 = mod(i,M)+1;
    dim2 = floor(i/M)+1;
    for k = 1:3
        y_dd_vec(i+1,k) = y_dd_block(dim1,dim2,k);
    end
end

% Finally, shuffle to make manual versions of y_tilde
y_tilde_cont = zeros(M*N,3);
for k = 1:3
    y_tilde_cont(:,k) = Gamma_MN' * y_dd_vec(:,k);
end

%% Direct Matrix Method
fprintf("Creating data using matrix method...\n");
% Create channel matrices
H = zeros(N*M,N*M,3);
for k = 1:3
    H(:,:,k) = gen_H_direct(T,N,M,Lp,Ln,chn_g,chn_tau,chn_v,ambig_vals(:,:,k),tap_t_indices,tap_f_indices);
end

% Create discrete versions of eta
eta_mat = zeros(N*M,3);
for k = 1:3
    eta_mat(:,k) = H(:,:,k) * sig_vec + w_vec(:,k);
end
eta_block = zeros(M,N,3);
for i = 0:N*M-1
    dim1 = mod(i,M)+1;
    dim2 = floor(i/M)+1;
    for k = 1:3
        eta_block(dim1,dim2,k) = eta_mat(i+1,k);
    end
end

% Convert to TF domain then DD
y_tf_block = zeros(M,N,3);
y_dd_block = zeros(M,N,3);
for k = 1:3
    y_tf_block(:,:,k) = F_M * eta_block(:,:,k);
    y_dd_block(:,:,k) = F_M' * y_tf_block(:,:,k) * F_N;
end

% Convert back to vectors
y_dd_vec = zeros(N*M,3);
for i = 0:N*M-1
    dim1 = mod(i,M)+1;
    dim2 = floor(i/M)+1;
    for k = 1:3
        y_dd_vec(i+1,k) = y_dd_block(dim1,dim2,k);
    end
end

% Finally, shuffle to make manual versions of y_tilde
y_tilde_mat = zeros(M*N,3);
for k = 1:3
    y_tilde_mat(:,k) = Gamma_MN' * y_dd_vec(:,k);
end

%% Comparison
fprintf("Comparing final results:\n")
% Directly compare samples at each filter
eta_diff = zeros(3,1);
y_tilde_diff = eta_diff;
for k = 1:3
    eta_diff(k) = norm(abs(eta_cont(:,k) - eta_mat(:,k)));
    y_tilde_diff(k) = norm(abs(y_tilde_cont(:,k) - y_tilde_mat(:,k)));
end

results = table(pulse_names,eta_diff,y_tilde_diff);
disp(results)
