clear; clc; close all;
%% Info
% Hi. 
% The best place to look for options on how to run this is to see the
% comments and constructors of each of the various classes. Options for
% some properties are also enumerations, so see Enums for more details there. 
% 
% Other helpful notes.
% * Everything is a class. 
% * Only "Signals" get passed between classes.
% * Every module must have a required_domain and required_fs property.
%   This is checked against all incoming data and then auto adjusted if
%   necessary. Even though some of theses can only be 1 possible value, I
%   made it so it has to be passed in, that way the experimenter is always
%   conscious of the DSP.
% 

%% Params
% Setup input signal
p.ofdm.n_scs = 1200;
p.ofdm.sc_spacing = 15e3; % 15kHz subcarrier spacing
p.ofdm.constellation = Constellation.QPSK;
p.ofdm.cp_length = 144; % Number of samples in cyclic prefix.
p.ofdm.n_symbols = 14;
p.ofdm.rms_input = 0.50; % digital RMS input.

% Setup DPD
p.dpd.postdistorter = 'GMP'; % Only option now. Though, hypothetically, another strategy could be added.
p.dpd.predistorter = 'GMP'; 
p.dpd.order = 9;
p.dpd.memory_depth = 4;
p.dpd.lag_depth = 0;  % 0 is a standard MP. >0 is GMP.
p.dpd.n_iterations = 8;
p.dpd.learning_rate = 0.5;
p.dpd.learning_method = LearningMethods.NEWTON; 
p.dpd.use_even = true; 
p.dpd.use_conj = 0;    % Conjugate branch. Currently only set up for MP (lag = 0)
p.dpd.use_dc_term = 0; % Adds an additional term for DC
p.dpd.required_domain = Domain.TIME;
p.dpd.required_fs = 200e6; % Needs to match the PA platform used below.

% Setup the PA simulator or TX board
PA_board = 'webRF'; % either 'WARP', 'webRF', or 'GMP'
p.pa.requried_fs = 200e6;
p.pa.required_domain = Domain.TIME;

%% Create the main modules
pa = PA.create(PA_board, p.pa);
dpd = ILA_DPD(p.dpd);

%% Run Expierement
tx_data = Signal.make_ofdm(p.ofdm);
w_out_dpd = board.transmit(tx_data.data);
dpd.perform_learning(tx_data.data, board);
w_dpd = board.transmit(dpd.predistort(tx_data.data));

%% Plot
w_out_dpd.plot_psd;
w_dpd.plot_psd;
dpd.plot_history;

switch PA_board
    case 'WARP'
        warp_params.n_boards = 1;        % Number of boards
        warp_params.RF_port  = 'A2B';    % Broadcast from RF A to RF B. Can also do 'B2A'
        requried_fs = 40e6;    % WARP board sampling rate.
    case 'none'
        board = PowerAmplifier(7, 4);
        requried_fs = 40e6;    % WARP board sampling rate.
    case 'webRF'
        dbm_power = -22;
        board = webRF(dbm_power);
end


% Setup OFDM

modulator = OFDM(ofdm);

% Create TX Data
[tx_data, ~] = modulator.use;
tx_data = Signal(tx_data, modulator.sampling_rate, rms_input);
tx_data.upsample(board.sample_rate)



