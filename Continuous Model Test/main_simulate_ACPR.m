%% MAIN_SIMULATE_ACPR
% This MATLAB file simulates doppler-delay domain information, produces a
% continuous signal passing through a rect/sinc/RRC filter, passes through
% a doubly-selective (t and f) channel, passes through another filter.
% The resulting sampled RX vector is then converted to the frequency domain
% and averaged across multiple iterations. The power spectral density is
% found, and from it the main band and out of band powers.
%
% Written by JRW, 5/20/2024

clear; clc; close all;
% System settings
res = 10;
q_vec = 20;
alpha = 1;
num_iters = 50;
Eb_N0_db = Inf;
M = 64;
N = 16;
filename_excel = sprintf('Saved Data\\OOBE_M%d_N%d_q%d-%d_alpha%.1f.xlsx',M,N,min(q_vec),max(q_vec),alpha);

% Set save spot
if ~exist("Saved Data", 'dir')
    mkdir("Saved Data")
end
pulse_names = ["Rectangular","Sinc",sprintf("RRC, α=%d",alpha)];
var_names = {'q','Main Band Power (dBW)', 'Adjascent Band Power (dBW)', 'ACPR (dB)'};
for k = 1:3
    writecell(var_names, filename_excel, 'Range', 'A1', 'Sheet', pulse_names(k));
end

% System initialization
Eb = 1;
N0 = Eb / (10^(Eb_N0_db / 10));
bit_order = [0,0;0,1;1,0;1,1];
S = Eb * [(sqrt(2)/2) + (1j*sqrt(2)/2);...
    (sqrt(2)/2) - (1j*sqrt(2)/2);...
    -(sqrt(2)/2) + (1j*sqrt(2)/2);...
    -(sqrt(2)/2) - (1j*sqrt(2)/2)];
T = 1/15000;
Ts = T/M;
f_nulls = [1,1/2,(1+alpha)/2] / Ts;
dt = Ts / res;
syms_per_f = M*N;
Fc = 4*10^9;
v = 500;

% Create needed matrices
Gamma_MN = gen_Gamma_MN(M,N);
F_M = gen_DFT(M);
F_N = gen_DFT(N);

for q_sel = 1:length(q_vec)
    q = q_vec(q_sel);
    filename_mat = sprintf('Saved Data\\OOBE_M%d_N%d_q%d_alpha%.1f.mat',M,N,q,alpha);
    filename_image = sprintf('Saved Figures\\PSD_M%d_N%d_q%d_alpha%.1f.png',M,N,q,alpha);

    % Sample from channel data to get channel length
    Lp =  q + floor((2510*10^(-9)) / Ts);
    Ln = -q;
    N0 = N0 * (N*M+Lp-Ln) / (N*M);

    % Time range setup
    t_range_sig = -Ts*(Lp):dt:N*T+Ts*(-Ln)-dt;
    t_range_filter = 0:dt:q*Ts;
    t_range_tx = -Ts*(Lp):dt:N*T+Ts*(-Ln+q)-dt;
    [~,tx_start_index] = min((t_range_tx+Ts*(Lp)).^2);
    [~,tx_end_index] = min((t_range_tx-(N*T+Ts*(-Ln)-dt)).^2);
    t_range_rx = -Ts*(Lp+q):dt:N*T+Ts*(-Ln+q)-dt;
    [~,rx_start_index] = min((t_range_rx+Ts*(Lp)).^2);
    [~,rx_end_index] = min((t_range_rx-(N*T+Ts*(-Ln)-dt)).^2);

    % Filter setup
    filters = zeros(3,length(t_range_filter));
    filters(1,:) = gen_rect_filter(t_range_filter,Ts,q/2);
    filters(2,:) = sinc_trunc((t_range_filter - q*Ts/2),Ts,q/2);
    filters(3,:) = RRCt_trunc((t_range_filter - q*Ts/2),alpha,Ts,q/2);

    % Normalize E of filters
    for k = 1:3
        filters(k,:) = filters(k,:) ./ sqrt(sum(filters(k,:).^2));
    end

    % Initialize needed matrices and vectors
    x_dd_block = zeros(M,N);
    sig_vec = zeros(N*M,1);
    sig_t = zeros((N*M+Lp-Ln)*res,1);
    s_t  = zeros(length(t_range_tx),3);
    r_t  = zeros(length(t_range_tx),3);
    nu_t = zeros(length(t_range_rx),3);
    w_t  = zeros(length(t_range_rx),3);

    % Define ACPR variables
    s_f = zeros(length(t_range_tx),3);
    r_f = zeros(length(t_range_tx),3);
    nu_f = zeros(length(t_range_rx),3);
    s_f_sum  = zeros(length(t_range_tx),3);
    r_f_sum  = zeros(length(t_range_tx),3);
    nu_f_sum = zeros(length(t_range_rx),3);

    fs2 = 1 / dt;
    N2 = (N*M-Ln+Lp+2*q)*res;
    bin_vals2 = 0:N2-1;
    N2_2 = ceil(N2/2);
    fax_Hz2 = (bin_vals2 - N2_2)*fs2/N2;

    % LOOP ------------------------------------------------------------------
    for iters = 1:num_iters
        % Data
        % Generate data
        [~,TX_2,x_dd_vec] = gen_data(bit_order,S,syms_per_f);
        TX_sym = Gamma_MN' * TX_2;
        x_tilde = Gamma_MN' * x_dd_vec;

        % Get relevent vectors and matrices from data
%         x_dd_block = zeros(M,N);
        for i = 0:N*M-1
            dim1 = mod(i,M)+1;
            dim2 = floor(i/M)+1;
            x_dd_block(dim1,dim2) = x_dd_vec(i+1);
        end
        x_tf_block = F_M * x_dd_block * F_N';
        sig_block = F_M' * x_tf_block;
%         sig_vec = zeros(N*M,1);
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
%         sig_t = zeros(numel(sig_buffer)*res,1);
        sig_t(1:res:end) = sig_buffer;

        % TX
        % Convolute time domain TX signal with first filter
%         s_t = zeros(length(t_range_tx),3);
        for k = 1:3
            s_t(:,k) = conv(filters(k,:),sig_t);
        end
        s_t_trunc = s_t(tx_start_index:tx_end_index,:);

        % Channel
        % Pass filtered signal through channel
        [chn_g,chn_tau,chn_v] = channel_generation(Fc,v);
        chn_tau = round(chn_tau/dt)*dt;

        r_t = zeros(length(t_range_tx),3);
        update_vals = floor(length(t_range_tx)*(linspace(.01,.99,34)));
        fprintf("Passing differently filtered signals through channel...\nProgress:                        |\n")
        for t_index = 1:length(t_range_tx)-1
            if ismember(t_index,update_vals)
                fprintf("x")
            end

            t = t_range_tx(t_index);

            rx_sums = zeros(3,1);
            for i = 1:length(chn_v)
                s_time = t - chn_tau(i);
                if s_time >= min(t_range_tx)
                    [~,s_index] = min(abs(t_range_tx - s_time).^2);
                else
                    s_index = length(t_range_tx);
                end
                for k = 1:3
                    sum_val = chn_g(i) * s_t(s_index,k) * exp(1j*2*pi*chn_v(i)*(t-chn_tau(i)));
                    rx_sums(k) = rx_sums(k) + sum_val;
                end
            end
            for k = 1:3
                r_t(t_index,k) = rx_sums(k);
            end
        end
        fprintf("\n");

        % Add noise to RX signal
        z_t = sqrt(N0/2) * (rand(length(r_t(:,1)),1) + 1j*rand(length(r_t(:,1)),1));
        r_t = r_t + z_t;
        r_t_trunc = r_t(tx_start_index:tx_end_index,:);

        % RX
        % Convolute time-domain RX signal and noise with second filter
%         nu_t = zeros(length(t_range_rx),3);
        w_t = nu_t;
        for k = 1:3
            nu_t(:,k) = conv(flip(filters(k,:)),r_t(:,k));
            w_t(:,k)  = conv(flip(filters(k,:)),z_t);
        end
        nu_t_trunc = nu_t(rx_start_index:rx_end_index,:);

        % Find power spectrum density of over-the-air signal
        for k = 1:3
            % Find frequency domain of current signal
            s_f(:,k) = abs(fftshift(fft(s_t(:,k))));
            r_f(:,k) = abs(fftshift(fft(r_t(:,k))));
            nu_f(:,k) = abs(fftshift(fft(nu_t(:,k))));
        end

        % Add frequency domain to sum for averaging
        s_f_sum = s_f_sum + s_f;
        r_f_sum = r_f_sum + r_f;
        nu_f_sum = nu_f_sum + nu_f;

        % Find PSD (square of average frequency response)
        psd_s = (s_f_sum / iters).^2;
        psd_r = (r_f_sum / iters).^2;
        psd_nu = (nu_f_sum / iters).^2;

        % % Plot PSD
        % figure(q_sel)
        % plot(fax_Hz2,10*log10(psd_nu));
        % title(sprintf("Power Spectral Density of Received Signals, Q=%d",q))
        % ylabel("Power/Frequency (dB/Hz)")
        % xlabel("Frequency (Hz)")
        % legend(pulse_names,Location="southwest")
        % grid on;
        % drawnow;
    end

    % % Plot PSD
    % figure(q_sel)
    % plot(fax_Hz2,10*log10(psd_nu));
    % title(sprintf("Power Spectral Density of Received Signals, Q=%d",q))
    % ylabel("Power/Frequency (dB/Hz)")
    % xlabel("Frequency (Hz)")
    % legend(pulse_names,Location="southwest")
    % grid on;

    % if ~exist("Saved Figures", 'dir')
    %     mkdir("Saved Figures")
    % end
    % saveas(gcf, filename_image);

    % Find power of each band given theoretical BW for each pulse
    df = fax_Hz2(2) - fax_Hz2(1);
    main_pwr = zeros(3,1);
    adj_pwr = zeros(3,1);
    out_pwr = zeros(3,1);
    for k = 1:3
        % f_null_sel = f_nulls(k);
        % [~,f_null_1p] = min(abs(fax_Hz2 - f_null_sel).^2);
        % [~,f_null_1n] = min(abs(fax_Hz2 + f_null_sel).^2);
        % [~,f_null_2p] = min(abs(fax_Hz2 - f_null_sel*2).^2);
        % [~,f_null_2n] = min(abs(fax_Hz2 + f_null_sel*2).^2);
        % 
        % main_pwr(k) = df * sum(psd_nu(f_null_1n:f_null_1p,k));
        % adj_pwr(k) = df * (sum(psd_nu(f_null_2n:f_null_1n-1,k)) + sum(psd_nu(f_null_1p+1:f_null_2p,k)));
        % out_pwr(k) = df * (sum(psd_nu(1:f_null_1n-1,k)) + sum(psd_nu(f_null_1p+1:end,k)));

        f_null_sel = f_nulls(k);
        [~,f0] = min(abs(fax_Hz2).^2);
        [~,f_null_1] = min(abs(fax_Hz2 - f_null_sel).^2);
        [~,f_null_2] = min(abs(fax_Hz2 - 3/Ts).^2);

        main_pwr(k) = df * sum(psd_nu(f0:f_null_1,k));
        adj_pwr(k) = df * sum(psd_nu(f_null_1:f_null_2,k));
    end

    % Calculate dB values
    inband_pwr_dB = 10*log10(main_pwr);
    adjband_pwr_dB = 10*log10(adj_pwr);
    % outband_pwr_dB = 10*log10(out_pwr);
    ACPR_dB = 10*log10(adj_pwr ./ main_pwr);
    % OOB_dB = 10*log10(out_pwr ./ main_pwr);

    % Make table of data and disply
    disp_data = [q*ones(3,1), inband_pwr_dB, adjband_pwr_dB, ACPR_dB];
    Table = array2table(disp_data);
    Table.Properties.RowNames = pulse_names;
    Table.Properties.VariableNames = var_names;
    disp(Table)

    % Write results to table
    for k = 1:3
        writecell(table2cell(Table(k,:)), filename_excel, 'Range', sprintf('A%d',q_sel+1), 'Sheet', pulse_names(k));
    end

    % Save source data for later
    save(filename_mat,"fax_Hz2","psd_nu");
end
