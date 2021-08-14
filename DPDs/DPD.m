classdef DPD < handle
    %DPD. Even though there is only 1 DPD in this repo now, I wanted to
    %implement this so that there would be a clear direction for scaling
    %up and a clear interface.
    
    properties
    end
    
    methods
        function obj = DPD()
        end
        
        function output_signal = learn(obj)
            data = input_signal.data;
            out_data = obj.subclass_learn(data);
            output_signal = Signal(out_data);
        end
        
        function output_signal = use(obj)
            data = input_signal.data;
            out_data = obj.subclass_use(data);
            output_signal = Signal(out_data);
        end
    end
    
    methods (Abstract)
        subclass_learn()
        subclass_use()
    end
    
    methods (Static)
        function obj = create(p)
            % Factory method for the DPD classes.
            subclass = str2func(p.subclass);
            obj = subclass(p);
            obj.required_domain = p.domain;
            obj.required_fs = p.required_fs;
        end
    end
end

