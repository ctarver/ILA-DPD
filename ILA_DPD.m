classdef ILA_DPD < handle
    %ILA_DPD. Inderect Learning Architecture DPD.
    %
    %  x(n)   +-----+ u(n) +-----+
    % +-----> | DPD +--+-> | PA  +------+
    %         +-----+  v   +-----+      |
    %                  --------+ e(n)   | y(n)
    %                  ^       |        |
    %                  |   +---v-+      |
    %                  +---+ DPD | <----+
    %               z(n)   +-----+
    %
    %
    %  MP DPD:
    %                +-----+    +---------------+
    %           +---->.    +---->b_1,1 ... b_1,M+-------+
    %           |    +-----+    +---------------+       |
    %           |                                       |
    %  x(n)     |    +-----+    +---------------+    +--v-+
    % +-------------->.|.|2+---->b_3,1 ... b_3,M+---->SUM +-------->
    %           |    +-----+    +---------------+    +--+-+
    %           |       .               .               ^
    %           |       .               .               |
    %           |       .                               |
    %           |    +-------+  +---------------+       |
    %           +---->.|.|p|1|  |b_p,1 ... b_p,M+-------+
    %                +-------+  +---------------+
    %
    %	Author:	Chance Tarver (2018)
    %		tarver.chance@gmail.com
    %
    
    properties
        order         % Nonlinear order of model. Can only be odd.
        use_even      % Include even order terms? true or false
        memory_depth  % Memory depth on each branch of the parallel hammerstein model
        lag_depth     % Memory depth of the lead/lag term in GMP
        nIterations   % Number of iterations used in the ILA learning
        coeffs        % DPD coefficients
        use_conj      % Use a conjugate branch as well
        use_dc_term   % use a dc term
        learning_rate % How much of the new iteration to use vs previous iteration. Should be in (0, 1]
        learning_method % Newton or ema
        coeff_history % Holds onto the coeffs used at each iteration
        result_history % Holds intermediate ACLR for each point during training in case of divergence.
    end
    
    methods
        function obj = ILA_DPD(params)
            %ILA_DPD. Make a DPD module
            
            if nargin == 0
                params.order = 7;
                params.memory_depth = 3;
                params.lag_depth = 0;
                params.nIterations = 3;
                params.use_conj = 0;
                params.use_dc_term = 0;
                params.learning_rate = 0.75;
            end
            
            if mod(params.order, 2) == 0
                error('Order of the DPD must be odd.');
            end
            
            obj.order = params.order;
            obj.use_even = params.use_even;
            obj.memory_depth = params.memory_depth;
            obj.lag_depth = params.lag_depth;
            obj.nIterations = params.nIterations;
            obj.learning_rate = params.learning_rate;
            obj.learning_method = params.learning_method;
            obj.use_conj = params.use_conj;
            obj.use_dc_term = params.use_dc_term;
            
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
        
        
        function perform_learning(obj, x, pa)
            %perform_learning. Perfrom ILA DPD.
            %
            % The PA output is the input to the postdistorter used for
            % learning. We want the error to be zero which happens when the
            % ouput of the pre and post distorters are equal. So we need:
            %
            %     e = 0
            % u - z = 0
            %     u = z
            %     u = Y * beta
            %
            % We can set this up as a least squares regression problem.
            
            obj.coeff_history = obj.coeffs;
            obj.result_history  = zeros(3, obj.nIterations+1);
            for iteration = 1:obj.nIterations
                % Forward through Predistorter
                u = obj.predistort(x);
                [y, test_signal] = pa.transmit(u); % Transmit the predistorted pa input
                obj.result_history(:, iteration) = test_signal.measure_all_powers;
                % Learn on postdistrter
                Y = setup_basis_matrix(obj, y);
                switch obj.learning_method
                    case 'newton'
                        post_distorter_out = obj.predistort(y);
                        error = u - post_distorter_out;
                        ls_result = ls_estimation(obj, Y, error);
                        obj.coeffs = obj.coeffs + (obj.learning_rate) * ls_result;
                    case 'ema'
                        ls_result = ls_estimation(obj, Y, u);
                        obj.coeffs = (1-obj.learning_rate) * obj.coeffs + (obj.learning_rate) * ls_result;
                end
                obj.coeff_history = [obj.coeff_history obj.coeffs];
            end
            % Need extra to evaluate final iteration
            u = obj.predistort(x);
            [~, test_signal] = pa.transmit(u); % Transmit the predistorted pa input
            obj.result_history(:, iteration+1) = test_signal.measure_all_powers;
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
        
        
        function number_of_coeffs = convert_order_to_number_of_coeffs(obj, order)
            %convert_order_to_number_of_coeffs. Helper function to easily
            %convert the order to number of coeffs. We need this because we
            %only model odd orders.
            
            if nargin == 1
                order = obj.order;
            end
            
            number_of_coeffs = (order + 1) / 2;
        end
        
        
        function out = predistort(obj, x)
            %predistort. Use the coeffs stored in object to predistort an
            %input.
            
            X = obj.setup_basis_matrix(x);
            out = X * obj.coeffs;
        end
        
        function plot_history(obj)
            % plot_history. Plots how the magnitude of the DPD coeffs
            % evolved over each iteration.
            figure(55);
            iterations = 0:obj.nIterations;
            subplot(1,2,1)
            hold on;
            plot(iterations, abs(obj.coeff_history'));
            title('History for DPD Coeffs Learning');
            xlabel('Iteration Number');
            ylabel('abs(coeffs)');
            grid on;
            subplot(1,2,2)
            
            plot(iterations, (obj.result_history'));
            grid on;
            title('Performance vs Iteration');
            ylabel('dBm');
            xlabel('Iteration Number');
            legend('L1', 'Main Channel', 'U1', 'Location', 'best')
        end
    end
end
