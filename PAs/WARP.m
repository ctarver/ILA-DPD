classdef WARP  < handle
    %WARP Class to control a WARPv3 board using WARPLab 7.7.1.
    %   This class wraps around the default interface for warplab so that it
    %   can be used in a more friendly way without having to worry about the
    %   proper syntax. It also helps to seperate all WARP related code from
    %   the rest of a project.
    
    properties
        nodes       % warp node objects
        node_tx     % transmitter board
        node_rx     % reciever board
        rx_length   % length of receive vectors
        tx_length
        ifc_ids
        eth_trig
        gains
        channels
        filters
        synchronization
        dc
        TX
        RX
    end
    
    methods
        function obj = WARP(params)
            %WARP class constructor. Sets up a WARP object with the warp board
            %configured with some gain on some channel.
            %
            % Args:
            %     N: Number of warpv3 nodes.
            %
            
            N = params.nBoards;
            
            USE_AGC = false;        % Use the AGC if running on WARP hardware
            
            switch params.RF_port
                case 'A2B'
                    obj.TX = 'RF_A';
                    obj.RX = 'RF_B';
                case 'B2A'
                    obj.TX = 'RF_B';
                    obj.RX = 'RF_A';
            end
            
            % RX variables
            obj.channels.RX        = 1;
            obj.gains.RxGainRF     = 1;  % Rx RF Gain in [1:3] (ignored if USE_AGC is true)
            obj.gains.RxGainBB     = 8; % Rx Baseband Gain in [0:31] (ignored if USE_AGC is true)
            obj.filters.RX_LPF     = 3;  % [0,1,2,3] for approx ![7.5,9.5,14,18]MHz corner
            obj.filters.RX_LPFFine = 2;  % Must be integer in [0,1,2,3,4,5] for [90,95,100,105,110]% scaling to LPF corner frequency
            obj.filters.RX_HPF     = 0;  % Must be 0 (HPF corner of 100 Hz) or 1 (default; HPF corner of 30 kHz) This filter setting is only used when RXHP is 'disable' (ie 0)
            obj.filters.RX_HP      = 'disable';  % (boolean MODE) MODE: true enables RXHP on the node when in manual gain control false disables RXHP on the node when in manual gain control
            
            % TX variables
            obj.channels.TX     = 1;
            obj.gains.TxGainBB  = 3;  % [0,1,2,3] for approx ![-5, -3, -1.5, 0]dB baseband gain
            obj.gains.TXGainRF  = 55; % [0:63] for approx [0:31]dB RF gain
            obj.filters.TX_LPF  = 2;  % [1,2,3] for approx [12,18,24]MHz corner frequencies ([24,36,48]MHz bandwidths)
            
            %DC offset
            obj.dc.real = 0;%0.0220;
            obj.dc.imag = 0;%-0.005;
            
            %Defaults till we compute
            obj.synchronization.delay = 40;
            obj.synchronization.phase = 0;
            obj.synchronization.done  = 0;
            obj.synchronization.sub_sample = 1;
            
            obj.nodes = wl_initNodes(N); %Setup N boards
            
            obj.node_tx = obj.nodes(1);  %Transmit on board 1
            %% PASS THIS OUT: obj.node_tx.serialNumber
            obj.node_rx = obj.nodes(1);  %Recieve  on board 1
            
            obj.eth_trig = wl_trigger_eth_udp_broadcast;
            wl_triggerManagerCmd(obj.nodes, 'add_ethernet_trigger', [obj.eth_trig]);
            trig_in_ids  = wl_getTriggerInputIDs(obj.nodes(1));
            trig_out_ids = wl_getTriggerOutputIDs(obj.nodes(1));
            wl_triggerManagerCmd(obj.nodes, 'output_config_input_selection', [trig_out_ids.BASEBAND], [trig_in_ids.ETH_A]);
            obj.ifc_ids = wl_getInterfaceIDs(obj.nodes(1));
            
            %Set Channel Settings
            wl_interfaceCmd(obj.nodes, obj.ifc_ids.(obj.TX), 'channel', 2.4, obj.channels.TX);
            wl_interfaceCmd(obj.nodes, obj.ifc_ids.(obj.RX), 'channel', 2.4, obj.channels.RX);
            
            %Set Filter Settings
            wl_interfaceCmd(obj.nodes, obj.ifc_ids.RF_ALL, 'rx_lpf_corn_freq', obj.filters.RX_LPF);
            wl_interfaceCmd(obj.nodes, obj.ifc_ids.RF_ALL, 'rx_lpf_corn_freq_fine', obj.filters.RX_LPFFine);
            wl_interfaceCmd(obj.nodes, obj.ifc_ids.RF_ALL, 'rx_hpf_corn_freq', obj.filters.RX_HPF);
            wl_interfaceCmd(obj.nodes, obj.ifc_ids.RF_ALL, 'rx_rxhp', obj.filters.RX_HP);
            wl_interfaceCmd(obj.nodes, obj.ifc_ids.RF_ALL, 'tx_lpf_corn_freq', obj.filters.TX_LPF);
            
            %Set Gain Settings
            wl_interfaceCmd(obj.nodes, obj.ifc_ids.RF_ALL, 'rx_gain_mode', 'manual');
            wl_interfaceCmd(obj.nodes, obj.ifc_ids.RF_ALL, 'rx_gains', ...
                obj.gains.RxGainRF, obj.gains.RxGainBB);
            wl_interfaceCmd(obj.nodes, obj.ifc_ids.RF_ALL, 'tx_gains', ...
                obj.gains.TxGainBB, obj.gains.TXGainRF);
            
            if(USE_AGC)
                wl_interfaceCmd(obj.node_rx, obj.ifc_ids.RF_ALL, 'rx_gain_mode', 'automatic');
                wl_basebandCmd(obj.nodes, 'agc_target', -13);
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Set up the Baseband parameters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Get the baseband sampling frequencies from the board
            ts_tx   = 1 / (wl_basebandCmd(obj.nodes(1), 'tx_buff_clk_freq'));
            ts_rx   = 1 / (wl_basebandCmd(obj.nodes(1), 'rx_buff_clk_freq'));
            ts_rssi = 1 / (wl_basebandCmd(obj.nodes(1), 'rx_rssi_clk_freq'));
            maximum_buffer_len = wl_basebandCmd(obj.node_tx, obj.ifc_ids.RF_A, 'tx_buff_max_num_samples');
            
            % Set the transmission / receptions lengths (in samples)
            %     See WARPLab user guide for maximum length supported by WARP hardware
            %     versions and different WARPLab versions.
            obj.tx_length    = 2^12;
            obj.rx_length    = obj.tx_length + 200;
            rssi_length  = obj.rx_length / (ts_rssi / ts_rx);
            
            % Check the transmission length
            if (obj.tx_length > maximum_buffer_len)
                error('Node supports max transmission length of %d samples.  Requested %d samples.', maximum_buffer_len, tx_length);
            end
            
            % Set the length for the transmit and receive buffers based on the transmission length
            wl_basebandCmd(obj.nodes, 'tx_length', obj.tx_length);
            wl_basebandCmd(obj.nodes, 'rx_length', obj.rx_length);
        end
        
        function out = transmit(obj,txData1,sync)
            %TRANSMIT Method for transmitting data on WARP board. Everything
            %should already be configured in the construction of the class
            %
            % Args:
            %     txData1: vector of data to send
            %     sync:    bool for performing sync or not.
            
            %Update the TX and RX lengths
            original_signal_input = txData1;
            obj.tx_length = length(txData1);
            obj.rx_length = obj.tx_length+100;
            wl_basebandCmd(obj.nodes, 'tx_length', obj.tx_length);
            wl_basebandCmd(obj.nodes, 'rx_length', obj.rx_length);
            
            %wl_basebandCmd(obj.nodes, 'continuous_tx', 1);
            if nargin == 2
                sync = 1;
            end
            
            %Add in DC Offset correction for WARP
            txData1 = txData1 - obj.dc.real - obj.dc.imag*j;
            
            %Normalize for TX
            max_real = max(abs(real(txData1)));
            max_imag = max(abs(imag(txData1)));
            max_max = max(max_real, max_imag);
            if max_max > 0.95
                warning('Saturating DAC of WARP.');
            end
            
            
            % Loop to get the right gain setting on RX.
            while(1)
                wl_basebandCmd(obj.node_tx, [obj.ifc_ids.(obj.TX)], 'write_IQ', txData1);
                
                %Enable TX/RX nodes and their buffers
                wl_interfaceCmd(obj.node_tx, obj.ifc_ids.(obj.TX), 'tx_en');
                wl_basebandCmd(obj.node_tx, obj.ifc_ids.(obj.TX), 'tx_buff_en');
                wl_interfaceCmd(obj.node_rx, obj.ifc_ids.(obj.RX), 'rx_en');
                wl_basebandCmd(obj.node_rx, obj.ifc_ids.(obj.RX), 'rx_buff_en');
                
                %Send trigger to start
                obj.eth_trig.send();
                
                %Get the values from the RX node
                rx_iq    = wl_basebandCmd(obj.node_rx, [obj.ifc_ids.(obj.RX)], 'read_IQ', 0, obj.rx_length);
                %rx_rssi  = wl_basebandCmd(obj.node_rx, [obj.ifc_ids.RF_B], 'read_RSSI', 0, rssi_length);
                
                if max(abs(real(rx_iq)))>0.99 | max(abs(imag(rx_iq)))>0.99
                    warning('Saturation in ADC. Reducing RX gain.');
                    obj.gains.RxGainBB = obj.gains.RxGainBB - 1;
                elseif max(abs(real(rx_iq))) < 0.5 | max(abs(imag(rx_iq))) < 0.5
                    warning('Underflow in ADC. Increaseing RX gain.');
                    obj.gains.RxGainBB = obj.gains.RxGainBB + 1;
                else
                    fprintf('The new RX Gain is %d.\n', obj.gains.RxGainBB);
                    break;
                end
                wl_interfaceCmd(obj.nodes, obj.ifc_ids.RF_ALL, 'rx_gains', obj.gains.RxGainRF, obj.gains.RxGainBB);
            end
            
            %Turn things off
            wl_basebandCmd(obj.nodes, obj.ifc_ids.(obj.TX), 'tx_rx_buff_dis');
            wl_interfaceCmd(obj.nodes, obj.ifc_ids.(obj.TX), 'tx_rx_dis');
            wl_basebandCmd(obj.nodes, obj.ifc_ids.(obj.RX), 'tx_rx_buff_dis');
            wl_interfaceCmd(obj.nodes, obj.ifc_ids.(obj.RX), 'tx_rx_dis');
            
            if(sync)
                if(obj.synchronization.done == 0) %Is this the 1st time?
                    % Perform synchronization.
                    obj = cyclosync(obj,rx_iq(obj.synchronization.delay:...
                        obj.synchronization.delay + length(txData1)-1),txData1);
                    obj.synchronization.done = 1 ; %Update flag
                    
                end
                rx_iq = rx_iq / obj.synchronization.phase;
                rx_iq = rx_iq(obj.synchronization.delay:obj.synchronization.delay + length(txData1)-1);
                
            end
            
            out = rx_iq * norm(original_signal_input) / norm(rx_iq);
            %out = rx_iq * sum(abs(original_signal_input)) / sum(abs(rx_iq));
            
            if  obj.synchronization.sub_sample
                %Set up a LS estimation for figuring out a subsample delay.
                X = [out [0; out(1:end-1)]];
                
                coeffs = (X'*X) \ (X'*original_signal_input);
                
                out = X*coeffs;
            end
            
        end
        function v = f_cfr_fitz(obj,z,L,N,Fs)
            %f_cfr_fitz Method for coarse CFO correction and estimation.
            
            R=zeros(1,N);
            
            for m=1:N,
                summa=0;
                for k=m+1:L,
                    summa=summa+z(k)*conj(z(k-m));
                end
                R(m)=1/(L-m)*summa;
            end
            
            v=1/(pi*N*(N+1)/Fs)*sum(atan2(imag(R),real(R)));
        end
        function obj = cyclosync(obj,x,y,varargin)
            %CYCLOSYNC	Synchronise two cyclic signals
            %	[DELAY,COEFF]=CYCLOSYNC(X,Y) finds the delay DELAY and complex
            %	coefficient COEFF (which accounts for gain and constant phase shift)
            %	to be applied to signal vector Y to synchronise with signal X. X and
            %	Y must be the same length, and are assumed to be one period of a
            %	periodic signal. In time domain, the corresponding continuous-time
            %	signal x(t) approximates COEFF * y(t-DELAY*Ts), where Ts is the
            %	sample interval.
            %
            %	[DELAY,COEFF,XSYNC]=CYCLOSYNC(X,Y) also returnes signal X
            %	synchronised to Y, i.e., delayed by -DELAY samples and divided by
            %	COEFF.
            %
            %	[DELAY,COEFF,YSYNC]=CYCLOSYNC(X,Y,...,'Y TO X') synchronises Y to X
            %	instead of X to Y by delaying Y by DELAY samples and multiplying
            %	by COEFF, and returns this in YSYNC.
            %
            %
            %	By delaying Y by DELAY samples and multiplying by COEFF, we obtain a
            %	signal synchronised with X. Or, by delaying X by -DELAY samples and
            %	dividing by COEFF, we obtain a signal synchronised with Y.
            %
            %
            %	Author:	Vesa Lehtinen (Aug 2012)
            %		ext-vesa.k.lehtinen@nokia.com
            %		vesa.lehtinen@tut.fi
            
            x=x(:);
            y=y(:);
            
            %----------.
            % Defaults |
            %----------'
            noScaling=true;
            sync='X TO Y';
            debug_corr=false;
            debug_phase=false;
            
            %---------.
            % Options |
            %---------'
            while numel(varargin)
                if ~ischar(varargin{1})
                    error 'Optional input args must start with an option name.'
                end
                switch upper(varargin{1})
                    case 'NO SCALING'
                        noScaling=true;
                        varargin(1)=[];
                    case {'X TO Y' 'Y TO X'}
                        sync=varargin{1};
                        varargin(1)=[];
                    case 'DEBUG:CORR'
                        debug_corr=true;
                        varargin(1)=[];
                    case 'DEBUG:PHASE'
                        debug_phase=true;
                        varargin(1)=[];
                    otherwise
                        error('Unknown option "%s".',varargin{1})
                end % switch
            end % while
            
            N=numel(x);
            if N ~= numel(y)
                error 'Signals must be the same length.'
            end
            
            f=[0:ceil(N/2)-1 -floor(N/2):-1]';
            
            X=fft(x);
            Y=fft(y);
            
            XY=X.*conj(Y);
            xya=abs(ifft(XY));
            delayInt=min(find(xya==max(xya)))-1;
            if debug_corr
                clf
                plot(0:N-1,xya)
                hold on
                plot([0 0]+delayInt,ylim,'m')
                hold off
                drawnow
                shg
                %%%%pause
            end % if debug_corr
            
            if 0
                %------------------------------.
                % Sync by interpolator fitting |
                %------------------------------'
                
                yd = y(1+mod((0:end-1)-delayInt, end));
                
                % yd ~ x + dx.*polyval(p1,n) + ddx.*polyval(p2,n) + dddx.*polyval(p3,n)
                % yd-x ~ [dx dx.*n... ddx ddx.*n... dddx dddx.*n...] * [p1 p2 p3]'
                
                n = (0:N-1)'-(N-1)/2;
                w = 2*pi*[0:ceil(N/2)-1 -floor(N/2):-1]'/N;
                pb = abs(w) <= 2*pi*1.02*28e6/81.6e6;
                dx = ifft(1j*w.*pb.*X);
                ddx = ifft(-w.^2.*pb.*X);
                dddx = ifft(-1j*w.^3.*pb.*X);
                
                if 0
                    I=2:N-1;
                else
                    I=1:N;
                end
                den = [dx n.*dx ddx n.*ddx n.^2.*ddx dddx n.*dddx n.^2.*dddx n.^3.*dddx];
                c = [real(den(I,:)); imag(den(I,:))] \ [real(yd(I)-x(I)); imag(yd(I)-x(I))]
                
                xBeforeFit = x;
                x = x + den * c;
                X = fft(x);
                
                if 01
                    cla
                    plot(abs(x-yd))
                    axis tight
                    ylim auto
                    title 'cyclosync: interpolator fitting'
                    if 0
                        beep
                        pause
                    else
                        drawnow
                    end
                end
                
                %error 'Sync by parameter fitting--UNIMPLEMENTED!'
                
            else
                xBeforeFit = [];
            end
            
            rot=@(n,x) 1+mod(n-1,numel(x));
            
            p=polyfit((-2:2)',xya(rot(1+delayInt+(-2:2),xya)),4);
            % Differentiate:
            pd=(4:-1:1).*p(1:4);
            % Find extrema:
            rt=roots(pd);
            % Remove invalid roots:
            rt(~~imag(rt) | rt<=-1 | rt>=1)=[];
            if isempty(rt)
                rt=0;
            end
            
            % Make sure there's only one root left, the one
            % that yields the highest polynomial value:
            pVal=polyval(p,rt);
            rt=rt(min(find(pVal==max(pVal))));
            
            % The total delay:
            delay=mod(delayInt+rt,N);
            closestDelay=delay-N*round(delay/N);
            
            if debug_corr
                hold on7
                plot([0 0]+delay,ylim,'r',delayInt+(-2:0.01:2),polyval(p,-2:0.01:2),'g')
                hold off
                drawnow
                shg
                pause
            end % if debug_corr
            
            
            switch upper(sync)
                
                case 'X TO Y'
                    synced=X.*exp(1j*2*pi*closestDelay*f/N);
                    if 01
                        fineDelay = -median(diff(angle(synced)-angle(Y)))/(2*pi);
                        delay=delay+fineDelay;
                        closestDelay=delay-N*round(delay/N);
                        if abs(fineDelay) > 0.1
                            warning 'Fine delay estimation failed.'
                            fineDelay=0;
                        else
                            synced=synced.*exp(1j*2*pi*fineDelay*f/N);
                        end % if
                    end % if 0[1]
                    coeff = Y \ synced;
                    if noScaling
                        coeff=coeff/abs(coeff);
                    end
                    synced=ifft(synced/coeff);
                    
                case 'Y TO X'
                    synced=Y.*exp(-1j*2*pi*closestDelay*f/N);
                    if 01
                        fineDelay = -median(diff(angle(X)-angle(synced)))/(2*pi);
                        delay=delay+fineDelay;
                        closestDelay=delay-N*round(delay/N);
                        if abs(fineDelay) > 0.1
                            warning 'Fine delay estimation failed.'
                            fineDelay=0;
                        else
                            synced=synced.*exp(-1j*2*pi*fineDelay*f/N);
                        end % if
                    end % if 0[1]
                    coeff = synced \ X;
                    if noScaling
                        coeff=coeff/abs(coeff);
                    end
                    synced=ifft(coeff*synced);
                    
                otherwise
                    error 'BUG!'
                    
            end % switch sync
            if debug_phase
                switch sync
                    case 'X TO Y'
                        plot(conj(Y).*synced)
                    case 'Y TO X'
                        plot(conj(synced).*X)
                end % switch
                axis tight
                axis equal
                hold on
                im=1j * 1e-6 * coeff * ~imag(coeff);  % Ensure complex plotting:
                plot([0 (coeff+im)/abs(coeff)*max(abs([xlim ylim]))],'r','LineWidth',3)
                hold off
                axis tight
                axis equal
                shg
                pause
            end % if debug_phase
            
            obj.synchronization.phase = coeff;
            obj.synchronization.delay = obj.synchronization.delay + round(delay);
            
        end
    end
end