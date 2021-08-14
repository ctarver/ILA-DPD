classdef ILA_DPD
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
        postdistorter   % GMP object used to perform actual predistortion.
        predistorter    % GMP object used to perform postdistortion to train the predistorter.
        n_iterations    % Number of iterations used in the ILA learning
        learning_rate   % How much of the new iteration to use vs previous iteration. Should be in (0, 1]
        learning_method % Newton or ema
        coeff_history   % Holds onto the coeffs used at each iteration
        result_history  % Holds intermediate ACLR for each point during training in case of divergence.
    end
    
    methods
        function obj = ILA_DPD(p)
            %ILA_DPD. Make a DPD module
            
            obj.n_iterations = p.n_iterations;
            obj.learning_rate = p.learning_rate;
            obj.learning_method = p.learning_method;
            
            obj.predistorter = GMP(p);
            obj.postdistorter = GMP(p);
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
        
        function plot_history(obj)
            % plot_history. Plots how the magnitude of the DPD coeffs
            % evolved over each iteration.
            figure(55);
            iterations = 0:obj.n_iterations;
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


classdef LearningMethods
   enumeration
       NEWTON, EMA
   end
end