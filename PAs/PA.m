classdef PA < handle
    
    properties
        
    end
    
    methods
        function obj = PA()
            
        end
        
        function output_signal = transmit(obj, input_signal)
           % Common interface for all subclasses. We extract the data, pass
           % to subclass, get result, wrap it in a signal, and return.
           % input_signal: Signal.
           
           data = input_signal.data;
           out_data = obj.subclass_transmit(data);
           output_signal = Signal(out_data);
        end
    end
    
    methods (Abstract)
        subclass_transmit(obj)
    end
    
    methods (Static)
        function obj = create(p)
            % Factory method for the PA classes.
            subclass = str2func(p.subclass);
            obj = subclass(p);
            obj.required_domain = p.required_domain;
            obj.required_fs = p.required_fs;
        end
    end 
end