classdef Signal < handle
    %Signal. This wrapper class is meant to encapsulate signals for easier
    %manipulation. This allows us to keep track of the sampling rate and to
    %allow us to easily up/down sample as needed as well as plot in the
    %time or frequency domain.
    %
    % Chance Tarver
    % tarver.chance@gmail.com
    % Feb 7, 2020
    
    properties
        data
        current_fs
        original_fs
        rms_power % dBM
        scale_factor
        resample_num
        resample_dem
        papr
        aclr
        obw % occupied bandwidth
    end
    
    methods
        function obj = Signal(data, current_fs, desired_rms)
            %Signal. Constructor for object. This will take in the signal,
            % its current sampling rate, and optionally a desired sampling
            % rms to rescale the signal to.
            % Inputs:
            %  - data: column vector of the signal.
            %  - current_fs: currect sampling rate of the data in Hz
            %  - desired_rms: optional. Will scale data to this rms.
            
            obj.current_fs = current_fs;
            obj.original_fs = current_fs;
            obj.data = data;
            if nargin == 3
                obj.normalize_to_this_rms(desired_rms);
            end
            obj.calculate_current_papr();
            obj.compute_occupied_bandwidth();
        end
        
        function normalize_to_this_rms(obj, this_rms)
            % Rescales the signal to a specific rms. Keeps track of the
            % scale factor in case you need to undo.
            obj.calculate_scale_factor(this_rms);
            obj.data = obj.data * obj.scale_factor;
            obj.calculate_current_rms_dbm();
            if abs(this_rms - obj.rms_power) > 0.01
                error('RMS is wrong.');
            end
        end
        
        function out = back_to_original_power(obj, in)
            % Undoes and scaling that has taken place.
            out = in / obj.scale_factor;
        end
        
        function set_new_resample_factors(obj, desired_rate)
            [obj.resample_num, obj.resample_dem] = rat(desired_rate/obj.current_fs);
        end
        
        function upsample(obj, desired_rate)
            % Upsamples to a desired sampling rate.
            if (nargin == 2)
                obj.set_new_resample_factors(desired_rate);
            end
            up = resample(obj.data, obj.resample_num, obj.resample_dem);
            obj.current_fs = obj.current_fs * obj.resample_num / obj.resample_dem;
            
            % Filter
            upsample_rate_inverse = obj.resample_dem/obj.resample_num;
            b = firls(50, [0 upsample_rate_inverse upsample_rate_inverse+0.1 1],[1 1 0 0]);
            up = [up; zeros(100,1)];
            obj.data = filter(b, 1, up);
        end
        
        function downsample(obj, desired_rate)
            if (nargin == 2)
                obj.set_new_resample_factors(desired_rate);
            end
            obj.data = resample(obj.data, obj.resample_dem, obj.resample_num);
            obj.current_fs = obj.current_fs * obj.resample_dem / obj.resample_num;
        end
        
        function append(obj, new_data)
            obj.data = [obj.data; new_data];
        end
        
        function calculate_current_rms_dbm(obj)
            obj.rms_power = 10*log10(norm(obj.data)^2/50/length(obj.data)) + 30;
        end
        
        function calculate_scale_factor(obj, desired_dbm_power)
            obj.scale_factor = sqrt(50 * length(obj.data) * 10 ^ ((desired_dbm_power-30)/10)) / norm(obj.data);
        end
        
        function calculate_current_papr(obj)
            obj.papr = 20*log10(max(abs(obj.data))*sqrt(length(obj.data))/norm(obj.data));
        end
        
        function plot_psd(obj, figure_handle)
            Nfft    = 1024;
            Window  = kaiser(1000, 9);
            Signal_PSD = 10*log10(fftshift(pwelch(obj.data, Window)));
            if nargin == 2
                plot(figure_handle, (-1:2/Nfft:1-2/Nfft)*((obj.current_fs)/(2e6)), Signal_PSD, ...
                    'LineWidth', 0.5);
            else
                figure(99);
                grid on;
                hold on;
                title('PSD');
                plot((-1:2/Nfft:1-2/Nfft)*((obj.current_fs)/(2e6)), Signal_PSD, ...
                    'LineWidth', 0.5);
                xlabel('Frequency (MHz)');
                ylabel('PSD (db/Hz)');
            end
        end
        
        function plot_time(obj)
            N = length(obj.data);
            period = 1/obj.current_fs;
            t = period * (1:N);
            plot(t, real(obj.data));
            xlabel('time (s)')
            ylabel('Magnitude');
        end
        
        function powers = measure_all_powers(obj)
            main_carrier_power = obj.measure_power('main')/50; % 50 ohm system.
            main_carrier_power_dbm = 10*log10(main_carrier_power/0.001);
            l1_power = obj.measure_power('L1')/50;
            l1_power_dbm = 10*log10(l1_power/0.001);
            u1_power = obj.measure_power('U1')/50;
            u1_power_dbm = 10*log10(u1_power/0.001);
            powers = [l1_power_dbm main_carrier_power_dbm u1_power_dbm]';
        end
        
        function powbp = measure_power(obj, type)
            if (obj.obw >= 15e6) && (obj.obw <= 20e6) % 20 MHz Signal
                ibw = 18e6;
                offset = 20e6;
            elseif (obj.obw >= 8e6) && (obj.obw <= 10e6) % 10 MHz Signal
                ibw = 9e6;
                offset = 10e6;
            elseif (obj.obw >= 3e6) && (obj.obw <= 5e6) % 5 MHz Signal
                ibw = 4500000;
                offset = 5e6;
            else
                warning('Unknown original fs in Signal.m measure_power method')
                powbp = -99;
                return;
            end
     
            switch type
                case 'main'
                    lower = -0.5*ibw;
                    upper = 0.5*ibw;
                case 'L1'
                    lower = -0.5*ibw - offset;
                    upper = 0.5*ibw - offset;
                case 'U1'
                    lower = -0.5*ibw + offset;
                    upper = 0.5*ibw + offset;
                otherwise
                    error('Unkown measure type in Signal.m, measure_power method');
            end
            F = [lower upper];
            
            try
                powbp = bandpower(obj.data, obj.current_fs, F);
            catch
                warning('Issue when measuring bandpower in Signal.m');
                powbp = -99;
            end
        end
        
        function compute_occupied_bandwidth(obj)
            obj.obw = obw(obj.data, obj.current_fs);
        end
    end
end

