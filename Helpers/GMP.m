classdef GMP < handle
    properties
        order         % Nonlinear order of model. Can only be odd.
        use_even      % Include even order terms? true or false
        memory_depth  % Memory depth on each branch of the parallel hammerstein model
        lag_depth     % Memory depth of the lead/lag term in GMP
        use_conj      % Use a conjugate branch as well
        use_dc_term   % use a dc term
        coeffs        % Model coefficients
    end
    
    methods
        function obj = GMP(p)
            obj.order = p.order;
            obj.use_even = p.use_even;
            obj.memory_depth = p.memory_depth;
            obj.lag_depth = p.lag_depth;
            obj.use_conj = p.use_conj;
            obj.use_dc_term = p.use_dc_term;
            
            if mod(p.order, 2) == 0
                error('Order of the DPD must be odd.');
            end
            
            % Start DPD coeffs being completely linear (no effect)
            if obj.use_even
                assert(obj.lag_depth == 0, 'GMP not yet supported for even terms. Set lag_depth=0');
                n_coeffs = obj.order * obj.memory_depth;
            else
                n_coeffs = obj.convert_order_to_number_of_coeffs * obj.memory_depth + ...
                    2*((obj.convert_order_to_number_of_coeffs-1) * obj.memory_depth * obj.lag_depth);
            end
            
            if obj.use_conj
                n_coeffs = 2*n_coeffs;
            end
            if obj.use_dc_term
                n_coeffs = n_coeffs + 1;
            end
            obj.coeffs = zeros(n_coeffs, 1);
            obj.coeffs(1) = 1;
        end

        function out = use(obj, x)
            %USE. Use the coeffs stored in object to apply the GMP to an
            %input.
            X = obj.setup_basis_matrix(x);
            out = X * obj.coeffs;
        end
        
        function number_of_coeffs = convert_order_to_number_of_coeffs(obj, order)
            %convert_order_to_number_of_coeffs. Helper function to easily
            %convert the order to number of coeffs. We need this because we
            %only model odd orders.
            
            if nargin == 1
                order = obj.order;
            end
            
            number_of_coeffs = (order + 1) / 2;
        end
        
        function beta = ls_estimation(obj, X, y)
            %ls_estimation
            % Solves problems where we want to minimize the error between a
            % lienar model and some input/output data.
            %
            %     min || y - X*beta ||^2
            %
            % A small regularlizer, lambda, is included to improve the
            % conditioning of the matrix.
            %
            
            % Trim X and y to get rid of 0s in X.
            X = X(obj.memory_depth+obj.lag_depth:end-obj.lag_depth, :);
            y = y(obj.memory_depth+obj.lag_depth:end-obj.lag_depth);
            
            lambda = 0.001;
            beta = (X'*X + lambda*eye(size((X'*X)))) \ (X'*y);
        end
        
        function X = setup_basis_matrix(obj, x)
            %setup_basis_matrix. Setup the basis matrix for the LS learning of
            %the PA parameters or for broadcasting through the PA model.
            %
            % obj.setup_basis_matrix(x)
            %
            % Inputs:
            %   x - column vector of the PA input signal.
            % Output:
            %   X - matrix where each column is the signal, delayed version of
            %   a signal, signal after going through a nonlinearity, or both.
            %
            %	Author:	Chance Tarver (2018)
            %		tarver.chance@gmail.com
            %
            
            number_of_basis_vectors = numel(obj.coeffs);
            X = zeros(length(x), number_of_basis_vectors);
            
            if obj.use_even
                step_size = 1;
            else
                step_size = 2;
            end
            
            % Main branch
            count = 1;
            for i = 1:step_size:obj.order
                branch = x .* abs(x).^(i-1);
                for j = 1:obj.memory_depth
                    delayed_version = zeros(size(branch));
                    delayed_version(j:end) = branch(1:end - j + 1);
                    X(:, count) = delayed_version;
                    count = count + 1;
                end
            end
            
            % Lag term
            for k = 3:step_size:obj.order  % Lag/Lead doesn't exist for k=1
                absolute_value_part_base = abs(x).^(k-1);
                for m = 1:obj.lag_depth
                    lagged_abs = [zeros(m,1); absolute_value_part_base(1:end-m)];
                    main_base = x .* lagged_abs;
                    for l = 1:obj.memory_depth
                        X(l:end, count) = main_base(1:(end-l+1));
                        count = count + 1;
                    end
                end
            end
            
            % Lead term
            for k = 3:step_size:obj.order  % Lag/Lead doesn't exist for k=1
                absolute_value_part_base = abs(x).^(k-1);
                for m = 1:obj.lag_depth
                    lead_abs = [absolute_value_part_base(1+m:end); zeros(m,1)];
                    main_base = x .* lead_abs;
                    for l = 1:obj.memory_depth
                        X(l:end, count) = main_base(1:(end-l+1));
                        count = count + 1;
                    end
                end
            end
            
            if obj.use_conj
                % Conjugate branch
                for i = 1:step_size:obj.order
                    branch = conj(x) .* abs(x).^(i-1);
                    for j = 1:obj.memory_depth
                        delayed_version = zeros(size(branch));
                        delayed_version(j:end) = branch(1:end - j + 1);
                        X(:, count) = delayed_version;
                        count = count + 1;
                    end
                end
            end
            
            % DC
            if obj.use_dc_term
                X(:, count) = 1;
            end
        end
    end
end