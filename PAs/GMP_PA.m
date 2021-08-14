classdef GMP_PA
    %GMP_PA Construct a PA, broadcast throough it, or model a PA.
    
    properties
        poly_coeffs  %matrix where each column corresponds to memory effects and each row is for a nonlinearity.
        order
        memory_depth
        nmse_of_fit
    end
    
    methods
        function obj = GMP_PA(order, memory_depth)
            %GMP_PA Construct an instance of this class.
            % This will initialize with a basic Parallel Hammerstein PA model
            % that was extracted from a WARP board.
            %
            % Args:
            %     order:        int with PA order. Should be odd. 1, 3, 5, etc.
            %     memory_depth: int with number of taps in FIR filter
            
            if nargin == 0
                order = 7;
                memory_depth = 4;
            end
            
            if mod(order,2) == 0
                error('Order must be odd.');
            end
            
            obj.order = order;
            obj.memory_depth = memory_depth;
            
            % Default polynomial coeffs derived from a WARP board.
            default_poly_coeffs = [ 0.9295 - 0.0001i, 0.2939 + 0.0005i, -0.1270 + 0.0034i, 0.0741 - 0.0018i;  % 1st order coeffs
                0.1419 - 0.0008i, -0.0735 + 0.0833i, -0.0535 + 0.0004i, 0.0908 - 0.0473i; % 3rd order coeffs
                0.0084 - 0.0569i, -0.4610 + 0.0274i, -0.3011 - 0.1403i, -0.0623 - 0.0269i;% 5th order coeffs
                0.1774 + 0.0265i, 0.0848 + 0.0613i, -0.0362 - 0.0307i, 0.0415 + 0.0429i]; % 7th order coeffs
            
            % Prune the model to have the desired number of nonlinearities and memory effects.
            obj.poly_coeffs = default_poly_coeffs(1:obj.convert_order_to_number_of_coeffs, 1:memory_depth);
        end
        
        
        function pa_output = transmit(obj, in)
            %transmit Broadcast the input data using the PA model currently stored in
            %the object.
            %
            %  obj.transmit(in) send in through the PH model that is stored in the
            %  object. It expands the input into a matrix where the columns
            %  are the different nonlinear branches or delayed versions of the
            %  nonlinear branches to model the FIR filter. A product can
            %  be done with the coefficient to get the PA output.
            %
            %	Author:	Chance Tarver (2018)
            %		tarver.chance@gmail.com
            %
            
            X = obj.setup_basis_matrix(in);
            coeffs = reshape(obj.poly_coeffs.',[],1);
            pa_output = X * coeffs;
        end
        
        
        function obj = make_pa_model(obj, in, out)
            %make_pa_model	Learn a PA model
            %	obj.make_pa_model(in, out) finds the best LS fit for a Parallel
            %	Hammerstein power amplifier. The 'in' is a column vector that is
            %	the signal that you put into the PA. The 'out' is a column
            %	vector that is the output of a real PA. This function will
            %	store the learned coefficients in obj.poly_coeffs. The PA model
            %	can be used by calling the transmit method.
            %
            %  This method also finds the NMSE of the derived PA model and
            %  stores it in obj.nmse_of_fit.
            %
            %  The LS regression solution is standard. Can be derrived by
            %  doing the sum_i [y_i - (beta_0 x_i + beta_! x_i)^2]
            %  optimization. The PA model is linear with respect to the
            %  coefficients.
            % 
            %  I am using a Regularization. It helps with the condition of the matrix
            %  http://www.signal.uu.se/Research/PCCWIP/Visbyrefs/Viberg_Visby04.pdf
            %  I just used a really small lambda.  
            %
            %	Author:	Chance Tarver (2018)
            %		tarver.chance@gmail.com
            %
            
            %% Construct signal matrix with basis vectors for each nonlinearity
            y = out;
            X = obj.setup_basis_matrix(in);
            
            %% LS solution to get the optimal coefficients.
            %coeffs = (X'*X) \ (X'*y); 
            lambda = 0.001;
            coeffs = (X'*X + lambda*eye(size((X'*X)))) \ (X'*y);
            
            %Reshape for easier to understand matrix of coeffs
            coeffs_transpose = reshape(coeffs, [obj.memory_depth, obj.convert_order_to_number_of_coeffs]);            
            obj.poly_coeffs = coeffs_transpose.';
            
            %% NMSE of the derived PA
            model_pa_output = obj.transmit(in);  
            obj.nmse_of_fit = obj.calculate_nmse(y, model_pa_output);
        end
        
        
        function nmse = calculate_nmse(~, desired, actual)
            %calculate_nmse. Calculate the normalized mean square error.
            % equivalent to sum (error)2 / sum(desired)^2
            nmse = norm(desired - actual)^2 / norm(desired)^2;
        end
        
        
        function X = setup_basis_matrix(obj, x)
            %setup_basis_matrix. Setup the basis matrix for the LS learning of
            %the PA parameters or for broadcasting through the PA model.
            %
            % obj.setup_basis_matrix(x)
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
    end
end
