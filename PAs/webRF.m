classdef webRF < PA
    %webRF Class wrapper for the webRF PA.
    % http://dpdcompetition.com/rfweblab/
    
    properties
        rms_in_dbm   % -24 seems to be a good value.
        rms_out_dbm
        Idc
        Vdc
        synchronization
        required_fs
        required_domain
    end
    
    methods
        function obj = webRF(p)
            %webRF Construct an instance of this class
            obj.rms_in_dbm = p.rms_in_dbm;
            
            assert(p.required_domain==Domain.TIME, 'WebRF Must be time domain');
            assert(p.required_fs==200e6, 'WebRF must 200 MHz sample rate');
        end
        
        function y = subclass_transmit(obj, x)
            %transmit. Take input signal, x, and broadcast it through the
            %RFWebLab PA.
            %
            %Args:
            %   -x: column vector. Will be normalized in RFWebLab function
            %
            %Returns:
            %   -y: column vector result from sending x through the PA. Y
            %       is normalized to be the same ||.||2 norm as x.
            
            if length(x) > 1000000
                warning("Too long for webRF.");
            end
            [y_original, obj.RMSout, obj.Idc, obj.Vdc] = RFWebLab_PA_meas_v1_1(x, obj.RMSin);
            
            % Need something to guarantee same as input length and aligned in TD.
            y = y_original(7:end);
            length_input = length(x);
            length_output = length(y);
            y = [y; zeros(length_input - length_output, 1)];
        end
    end
end