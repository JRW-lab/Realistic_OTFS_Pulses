function [TX_bit,TX_sym,s] = gen_data(bin_order,S,num_data_syms)
% Generate random bits, symbol equivalents, and modulated symbols

% Create TX symbols
TX_sym = randi(length(S),num_data_syms,1)-1;

% Create TX bits
TX_bit = bin_order(TX_sym+1,:);

% Create modulated symbols s
s = S(TX_sym+1);