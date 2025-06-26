
clear; clc; close all;
M = 64;
N = 16;
q = 8;
alpha = 1;
pulse_names = ["Rectangular","Sinc",sprintf("RRC, α=%d",alpha)];
T = 1/15000;
Ts = T/M;
legend_vec = ["r-","k--","b-."];

fig_num = 1;
font_val = 16;
line_val = 2;

filename_mat = sprintf('Saved Data\\OOBE_M%d_N%d_q%d_alpha%.1f.mat',M,N,q,alpha);
loaded_file = load(filename_mat);
fax_Hz2 = loaded_file.fax_Hz2;
psd_nu = loaded_file.psd_nu;

[~,low_f]  = min(abs(fax_Hz2 + 3/Ts).^2);
[~,high_f] = min(abs(fax_Hz2 - 3/Ts).^2);

figure(fig_num)
plot(fax_Hz2(low_f:high_f)*Ts,10*log10(psd_nu(low_f:high_f,1)),legend_vec(1));
hold on
plot(fax_Hz2(low_f:high_f)*Ts,10*log10(psd_nu(low_f:high_f,2)),legend_vec(2));
plot(fax_Hz2(low_f:high_f)*Ts,10*log10(psd_nu(low_f:high_f,3)),legend_vec(3));
% title(sprintf("Power Spectral Density of Received Signals, Q=%d",q))
ylabel("Power/Frequency (dB/Hz)")
xlabel("f (in 1/T_s)")
legend(pulse_names,Location="southwest")
grid on;


% Set font size
set(gca, 'FontSize', font_val);

% Adjust plotted line and border line size
lines = findall(gcf, 'Type', 'Line');
set(lines, 'LineWidth', line_val);
ax = gca;
ax.LineWidth = line_val;

% Save figure
saveas(figure(fig_num), sprintf('Figure%d_%s.eps',fig_num,string(datetime('now', 'format', 'MM.dd.uuuu_HH.mm'))), 'epsc');