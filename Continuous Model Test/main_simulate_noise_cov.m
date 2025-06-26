%% MAIN_SIMULATE_NOISE_COV
% This MATLAB file simulates AWGN and passes it through filters before
% generating the covariance matrix to determine if the result is still
% white.
%
% Written by JRW, 5/20/2024

clear; clc; close all;
if ~exist("Saved Data", 'dir')
    mkdir("Saved Data")
end

% System settings
res_discretization = 10;
q = 2;
alpha = 1;
Eb_N0_db = 0;
M = 64;
N = 16;
num_trials = 1000;
tap_range = 10;

% System initialization
Eb = 1;
N0 = (Eb / (10^(Eb_N0_db / 10)));
T = 1/15000;
Ts = T/M;
dt = Ts / res_discretization;
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

% Filters
% Filter setup
filters = zeros(3,length(t_range_filter));
filters(1,:) = gen_rect_filter(t_range_filter,Ts,q/2);
filters(2,:) = sinc_trunc((t_range_filter - q*Ts/2),Ts,q/2);
filters(3,:) = RRCt_trunc((t_range_filter - q*Ts/2),alpha,Ts,q/2);

% Normalize E of filters
for k = 1:3
    filters(k,:) = filters(k,:) ./ sqrt(sum(filters(k,:).^2));
end

% Run trials to simulate noise covariance
update_vals = 0:0.01:1;
fprintf("Passing differently filtered signals through channel...\nProgress:                                                                               |\n")
w_samples = zeros(N*M,num_trials,3);
z_tilde_samples = zeros(N*M,num_trials,3);
R_w_sum = zeros(N*M,N*M,3);
R_w = zeros(N*M,N*M,3);
R_zddt_sum = zeros(N*M,N*M,3);
R_zddt = zeros(N*M,N*M,3);
for g = 1:num_trials
    if ismember(g/num_trials,update_vals)
        fprintf("x")
    end

    % Generate noise
    z_t = sqrt(N0/2) * (randn(length(t_range_tx),3) + 1j*randn(length(t_range_tx),3));
    % z_t = sqrt(N0/2) * [zeros(res_discretization*Lp+10,3); ones(1,3); zeros(length(t_range_tx)-res_discretization*Lp-11,3)];

    % Convolute time-domain RX signal and noise with second filter
    w_t = zeros(length(t_range_rx),3);
    for k = 1:3
        w_t(:,k) = conv(flip(filters(k,:)),z_t(:,k));
    end

    % Sample RX signal after filtering for discrete representation
    w_vec_ext = zeros(N*M+Lp-Ln+2*q,3);
    for k = 1:3
        w_vec_ext(:,k) = w_t(1:res_discretization:end,k);
    end

    % Truncate sampled vector to just one frame
    w_vec = zeros(N*M,3);
    for k = 1:3
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
    z_tf_vec = zeros(N*M,3);
    z_tilde = zeros(N*M,3);
    for k = 1:3
        z_tf_block = F_M * w_block(:,:,k);
        for i = 0:N*M-1
            dim1 = mod(i,M)+1;
            dim2 = floor(i/M)+1;
            z_tf_vec(i+1,k) = z_tf_block(dim1,dim2);
        end
        z_dd_block = F_M' * z_tf_block * F_N;
        z_dd_vec = zeros(N*M,1);
        for i = 0:N*M-1
            dim1 = mod(i,M)+1;
            dim2 = floor(i/M)+1;
            z_dd_vec(i+1) = z_dd_block(dim1,dim2);
        end
        z_tilde(:,k) = Gamma_MN' * z_dd_vec;
    end

    % Sample the new noise vectors
    for k = 1:3
        w_samples(:,g,k) = w_vec(:,k);
        z_tilde_samples(:,g,k) = z_tilde(:,k);

        R_w_sum(:,:,k) = R_w_sum(:,:,k) + w_vec(:,k) * w_vec(:,k)';
        R_zddt_sum(:,:,k) = R_zddt_sum(:,:,k) + z_tilde(:,k) * z_tilde(:,k)';
        if g > 1
            R_w(:,:,k) = R_w_sum(:,:,k) / ((g-1)*N0);
            R_zddt(:,:,k) = R_zddt_sum(:,:,k) / ((g-1)*N0);

            % figure(1)
            % subplot(3,1,k)
            % plot(abs(R_zddt(1,:,k)));
            % drawnow
        end
    end
end
fprintf("\nComplete!\n");

R_w(1,1:3,2)

%%

% Create average covariance diagonal values
R_sum = zeros(N*M,N*M,3);
% R_sum2 = zeros(N*M,N*M,3);
for k = 1:3
    R_sum(:,:,k) = z_tilde_samples(:,:,k) * z_tilde_samples(:,:,k)';
    R_sum2(:,:,k) = w_samples(:,:,k) * w_samples(:,:,k)';
end
R = (R_sum / (num_trials-1)) / N0;

R2 = (R_sum2 / (num_trials-1)) / N0;
% for k = 1:3
%     R2(:,:,k) = kron(eye(N),F_M) * R2(:,:,k) * kron(eye(N),F_M)';
%     R2(:,:,k) = kron(F_N',F_M)' * R2(:,:,k) * kron(F_N',F_M);
%     R2(:,:,k) = Gamma_MN' * R2(:,:,k) * Gamma_MN;
% end

% Create average covariance diagonal values
cov_vals = zeros(N*M,3);
for k = 1:3
    for l = 0:tap_range-1
        cov_vals(l+1,k) = mean(diag(R(:,:,k),l));
    end
    % cov_vals(:,k) = cov_vals(:,k) / cov_vals(1,k);
end

% Generate covariance
RzDD_tilde = zeros(N*M,N*M,3);
D = zeros(N*M,N*M,3);
for k = 1:3
    RzDD_tilde(:,:,k) = toeplitz(cov_vals(:,k)) / N0;

    [V_z,Omega_z,~] = svd(RzDD_tilde(:,:,k));
    D(:,:,k) = sqrt(inv(Omega_z)) * V_z';

    % figure(k)
    % mesh(abs(RzDD_tilde(:,:,k)))
end

filename_mat = sprintf('Saved Data\\NWM_N%d_M%d_alpha%.f_q%d.mat',M,N,alpha,q);
save(filename_mat,"RzDD_tilde","D");
