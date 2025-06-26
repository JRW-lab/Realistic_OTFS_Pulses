function H = gen_H_direct(T,N,M,Lp,Ln,chn_g,~,chn_v,ambig_table,tap_t_range,tap_f_range)
% This generates the time-domain channel matrix for an OTFS system
% INPUTS:
%   Fc: carrier frequency
%   v:  maximum vehicle velocity
%   df: subcarrier spacing
%   N:  number of time symbols
%   M:  number of subcarriers
%
% Coded by Jeremiah Rhys Wimer, 3/24/2024

% rng('default');

% % Debugging inputs
% clear; clc;
% obj = comms_obj_OTFS;
% Fc = obj.Fc;
% v = obj.v_vel;
% df = obj.sbcar_spacing;
% N = 4;
% M = 4;
% shape = "rrc";
% alpha = 1;
% q = 1;

% Internal settings
Ts = T / M;

% Check to make sure channel length is still valid
if M <= min(Lp - Ln + 1)
    error("M MUST BE GREATER THAN TOTAL CHANNEL LENGTH")
end

% % METHOD 1 - USES FOR LOOPS
% % Create channel coefficients
% h = zeros(M*N,Lp-Ln+1);
% for k = 0:(N*M)-1
%     count = 0;
%     for l = Ln:Lp
%         count = count + 1;
%         % Render a vector of exponential values
%         exp_vals = exp(1j*2*pi*chn_v*(k-l)*Ts);
% 
%         % Render a vector of cross-ambiguity function values
%         tap_t_range = l*Ts - chn_tau;
%         tap_f_range = chn_v;
% 
%         % Normalize t and f ranges by dt and df, making sure all entries are valid
%         % indices
%         dt = ambig_t_range(2) - ambig_t_range(1);
%         df = ambig_f_range(2) - ambig_f_range(1);
%         tap_t_range = round(tap_t_range ./ dt) + ceil(length(ambig_table(:,1))/2);
%         tap_t_range(tap_t_range < 1) = 1;
%         tap_t_range(tap_t_range > length(ambig_table(:,1))) = length(ambig_table(:,1));
%         tap_f_range = round(tap_f_range ./ df) + ceil(length(ambig_table(1,:))/2);
%         tap_f_range(tap_f_range < 1) = 1;
%         tap_f_range(tap_f_range > length(ambig_table(1,:))) = length(ambig_table(1,:));
%         % Gen ambiguity values
%         ambig_vals = tap_t_range * 0;
%         for m = 1:length(tap_t_range)
%             try
%                 ambig_vals(m) = ambig_table(tap_t_range(m), tap_f_range(m));
%             catch
%                 1;
%             end
%         end
% 
% %         ambig_vals = 0 * chn_g;
% %         for p = 1:length(tap_t_range)
% %             ambig_vals(p) = ambig_direct(tap_t_range(p),tap_f_range(p),Ts,shape,alpha,q,res_ambig);
% %         end
% 
%         % Create sum vector
%         sum_vals = chn_g .* exp_vals .* ambig_vals;
% 
%         % Sum all elements to make value of h
%         h(k+1,count) = sum(sum_vals,"all");
%     end
% end

% % METHOD 2 - USES 2D SPACE
% % Define vector of l given L- and L+
% l = (Ln:Lp).';
% 
% % Premake all possible t and f ranges for h
% ambig_t_range = l*Ts - chn_tau;
% ambig_f_range = ones(Lp-Ln+1,1) .* chn_v;
% 
% % Create channel coefficients
% h = zeros(M*N,Lp-Ln+1);
% for k = 0:(N*M)-1
%     % Render a vector of exponential values
%     exp_vals = exp(1j.*2.*pi.*chn_v.*(k-l).*Ts);
% 
%     % Render a vector of cross-ambiguity function values
%     ambig_vals = ambig(ambig_t_range, ambig_f_range, Ts, shape, alpha);
% 
%     % Create sum vector
%     sum_vals = chn_g .* exp_vals .* ambig_vals;
% 
%     % Sum all elements to make value of h
%     h(k+1,:) = sum(sum_vals,2).';
% end

% METHOD 3 - USES 3D SPACE
% Define vector of l given L- and L+
l = (Ln:Lp).';

% % Make all possible t and f ranges for h
% tap_t_range = (l*Ts - chn_tau) .* ones(Lp-Ln+1,length(chn_g),N*M);
% ambig_f_range = (ones(Lp-Ln+1,1) .* chn_v) .* ones(Lp-Ln+1,length(chn_g),N*M);
% 
% % Normalize t and f ranges by dt and df, making sure all entries are valid
% % indices
% tap_t_range = round(tap_t_range ./ dt) + ceil(length(ambig_table(:,1))/2); 
% tap_t_range(tap_t_range < 1) = 1;
% tap_t_range(tap_t_range > length(ambig_table(:,1))) = length(ambig_table(:,1));
% ambig_f_range = round(ambig_f_range ./ df) + ceil(length(ambig_table(1,:))/2); 
% % ambig_f_range(ambig_f_range < 1) = 1;
% % ambig_f_range(ambig_f_range > length(ambig_table(1,:))) = length(ambig_table(1,:));

% Make all possible exponential values for h
k = reshape(0:(N*M-1), [1, 1, N*M]);
exp_vals = exp(1j.*2.*pi.*chn_v.*(k-l).*Ts);

% Make ambiguity values and sum to make h
% ambig_vals = ambig(ambig_t_range, ambig_f_range, Ts, shape, alpha);
% ambig_vals = ambig_direct(ambig_t_range, ambig_f_range, Ts, shape, alpha, q);
ambig_vals = 0 * tap_t_range;
for k = 1:size(tap_t_range,1)
    for l = 1:size(tap_t_range,2)
        for m = 1:size(tap_t_range,3)
            try
                ambig_vals(k,l,m) = ambig_table(tap_t_range(k,l,m),tap_f_range(k,l,m));
            catch
                error("AMBIG TABLE ERROR")
            end
        end
    end
end
sum_vals = chn_g .* exp_vals .* ambig_vals;
h = squeeze(sum(sum_vals,2)).';

% Create channel matrix from coefficients
H = zeros(M*N,M*N,1);
H(:,1:Lp-Ln+1) = fliplr(h(:,1:Lp-Ln+1));
for k = 1:M*N
    H(k,:) = circshift(H(k,:),k-Lp-1);
end
