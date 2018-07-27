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
    %  PH DPD:
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
        memory_depth  % Memory depth on each branch of the parallel hammerstein model
        nIterations   % Number of iterations used in the ILA learning
        block_size    % Block size used for each iteration in the learning
        coeffs        % DPD coefficients
    end
    
    methods
        function obj = ILA_DPD(params)
            %ILA_DPD. Make a DPD module
            
            if nargin == 0
                params.order = 7;
                params.memory_depth = 3;
                params.nIterations = 3;
                params.block_size = 50000;
            end
            
            if mod(params.order, 2) == 0
                error('Order of the DPD must be odd.');
            end
            
            obj.order = params.order;
            obj.memory_depth = params.memory_depth;
            obj.nIterations = params.nIterations;
            obj.block_size = params.block_size;
            
            % Start DPD coeffs being completely linear (no effect)
            obj.coeffs = zeros(obj.convert_order_to_number_of_coeffs, obj.memory_depth);
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
            
            for iteration = 1:obj.nIterations
                % Forward through Predistorter
                u = obj.predistort(x);
                y = pa.transmit(u); % Transmit the predistorted pa input
                
                % Learn on postdistrter
                Y = setup_basis_matrix(obj, y);
                beta = ls_estimation(obj, Y, u);
                obj.coeffs = vector_to_matrix_coeffs(obj, beta);
            end
        end
        
        
        function out = vector_to_matrix_coeffs(obj, in)
            %Reshape for easier to understand matrix of coeffs. Input the
            %vector form for LS formulation, and the function will output a
            %reshaped version where each row corresponds to a different
            %nonlinear branch and each column is a different memory depth.
            
            coeffs_transpose = reshape(in, [obj.memory_depth, obj.convert_order_to_number_of_coeffs]);
            out = coeffs_transpose.';
        end
        
        
        function beta = ls_estimation(~, X, y)
            %ls_estimation
            % Solves problems where we want to minimize the error between a
            % lienar model and some input/output data.
            %
            %     min || y - X*beta ||^2
            %
            % A small regularlizer, lambda, is included to improve the
            % conditioning of the matrix.
            %
            
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
            
            number_of_basis_vectors = obj.memory_depth * obj.convert_order_to_number_of_coeffs;
            X = zeros(length(x), number_of_basis_vectors);
            
            count = 1;
            for i = 1:2:obj.order
                branch = x .* abs(x).^(i-1);
                for j = 1:obj.memory_depth
                    delayed_version = zeros(size(branch));
                    delayed_version(j:end) = branch(1:end - j + 1);
                    X(:, count) = delayed_version;
                    count = count + 1;
                end
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
            beta = reshape(obj.coeffs.', [], 1);
            out = X * beta;
        end
    end
end

