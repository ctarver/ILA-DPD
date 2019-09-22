clear; clc; close all;
%% Setup Everything

% Add the submodules to path
addpath(genpath('OFDM-Matlab'))
addpath(genpath('WARPLab-Matlab-Wrapper'))
addpath(genpath('Power-Amplifier-Model'))

rms_input = 0.50;

% Setup the PA simulator or TX board
PA_board = 'webRF'; % either 'WARP', 'webRF', or 'none'
switch PA_board
    case 'WARP'
        warp_params.nBoards = 1;         % Number of boards
        warp_params.RF_port  = 'A2B';    % Broadcast from RF A to RF B. Can also do 'B2A'
        board = WARP(warp_params);
        Fs = 40e6;    % WARP board sampling rate.
    case 'none'
        board = PowerAmplifier(7, 4);
        Fs = 40e6;    % WARP board sampling rate.
    case 'webRF'
        dbm_power = -26;
        board = webRF(dbm_power);
        Fs = 200e6;   % webRF sampling rate.
end

% Setup OFDM
ofdm_params.nSubcarriers = 600;
ofdm_params.subcarrier_spacing = 15e3; % 15kHz subcarrier spacing
ofdm_params.constellation = 'QPSK';
ofdm_params.cp_length = 144; % Number of samples in cyclic prefix.
ofdm_params.nSymbols = 10;
modulator = OFDM(ofdm_params);

% Create TX Data
[tx_data, ~] = modulator.use;
upsampled_tx_data = up_sample(tx_data, Fs, modulator.sampling_rate);
tx_data = normalize_for_pa(upsampled_tx_data, rms_input);

% Setup DPD
dpd_params.order = 11;
dpd_params.memory_depth = 4;
dpd_params.nIterations = 2;
dpd_params.block_size = 50000;

dpd_params.use_conj = 1;
dpd_params.use_dc_term = 1;
conj_dpd = ILA_DPD(dpd_params);

dpd_params.use_conj = 0;
dpd_params.use_dc_term = 0;
dpd = ILA_DPD(dpd_params);

%% Run Expierement
w_out_dpd = board.transmit(tx_data);
conj_dpd.perform_learning(tx_data, board);
w_conj_dpd = board.transmit(conj_dpd.predistort(tx_data));

dpd.perform_learning(tx_data, board);
w_dpd = board.transmit(dpd.predistort(tx_data));

%% Plot
plot_results('psd', 'Original TX signal', tx_data, Fs)
plot_results('psd', 'No DPD', w_out_dpd, Fs)
plot_results('psd', 'With Normal DPD', w_dpd, Fs)
plot_results('psd', 'With Conjug DPD', w_conj_dpd, Fs)

%% Some helper functions
function out = up_sample(in, Fs, sampling_rate)
upsample_rate = floor(Fs/sampling_rate);
up = upsample(in, upsample_rate);
b = firls(255,[0 (1/upsample_rate -0.02) (1/upsample_rate +0.02) 1],[1 1 0 0]);
out = filter(b,1,up);
%beta = 0.25;
%upsample_span = 60;
%sps = upsample_rate;
%upsample_rrcFilter = rcosdesign(beta, upsample_span, sps);
%out = upfirdn(in, upsample_rrcFilter, upsample_rate);
end

function [out, scale_factor] = normalize_for_pa(in, RMS_power)
scale_factor = RMS_power/rms(in);
out = in * scale_factor;
if abs(rms(out) - RMS_power) > 0.01
    error('RMS is wrong.');
end

max_real = max(abs(real(out)));
max_imag = max(abs(imag(out)));
max_max = max(max_real, max_imag);
fprintf('Maximum value: %1.2f\n', max_max);
end
