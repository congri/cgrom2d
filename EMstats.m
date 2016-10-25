classdef EMstats
    %Class describing the Expectation-Maximization procedure and storing the data sequentially
    
    properties
        %% Optimization data
        MCMCStepWidth
        
        %% Optimal params, theta_cf
        W
        S
        mu
        %% Optimal params, theta_c
        theta
        sigma
        
    end
    
    properties (SetAccess = private)
        %read-only properties
        
        %Optimization params
        maxIterations
        
    end
    
    methods
        
        function obj = setMaxIterations(obj, maxIterations)
            
            %% set value
            obj.maxIterations = int32(maxIterations);
        end
        
        function obj = prealloc(obj, fineData, domainf, domainc, nBasis)
            %preallocation of data arrays
            
            %% for data record, single precision is sufficient
            obj.MCMCStepWidth = zeros(fineData.nSamples, obj.maxIterations + 1, 'single');
            w{1} = spalloc(domainf.nNodes, domainc.nNodes, 4*domainf.nNodes);
            obj.W = repmat(w, obj.maxIterations + 1, 1);
%             obj.W = zeros(domainf.nNodes, domainc.nNodes, obj.maxIterations, 'single');
            obj.theta = zeros(nBasis, obj.maxIterations + 1, 'single');
            obj.sigma = zeros(obj.maxIterations + 1, 1, 'single');
            s{1} = ones(domainf.nNodes, 1);
            obj.S = repmat(s, obj.maxIterations + 1, 1);
            obj.mu = zeros(domainf.nNodes, obj.maxIterations, 'single');
        end
        
        function plotData(obj)
            %further plotting functions
            addpath(genpath('~/matlab/toolboxes/plotting'))
            %plot the final data
            f = figure;
            set(f, 'Position', [680 678 900 700]) %[position x, psition y, width, height]
            subplot(2, 2, 1)
            iterations = 1:(obj.maxIterations + 1);
            p_theta = plot(iterations, obj.theta);
            set(gca, 'fontsize', 15)
            set(p_theta, 'linewidth', 2)
            %compute perceptually distinct colors
            c_theta = distinguishable_colors(size(obj.theta, 1));
            for i = 1:numel(p_theta)
                set(p_theta(i), 'color', c_theta(i, :))
            end
            xlabel('iteration i')
            ylabel('\theta_c')
            xlim([1 obj.maxIterations + 1])
            axis square
            grid on
%             ax = gca;
%             ax.XMinorGrid = 'on';
%             ax.YMinorGrid = 'on';
            
            
            subplot(2, 2, 2)
            p_sigma = semilogy(iterations, obj.sigma);
            set(gca, 'fontsize', 15)
            set(p_sigma, 'linewidth', 2)
            c_sigma = distinguishable_colors(1);
            set(p_sigma, 'color', c_sigma)
            axis square
            grid on
            xlabel('iteration i')
            ylabel('\sigma')
            xlim([1 obj.maxIterations + 1])
            ax = gca;
            ax.YMinorTick = 'on';
            ax.XMinorGrid = 'on';
            ax.YMinorGrid = 'on';
            
            
            subplot(2, 2, 3)
            Splt = zeros(size(obj.S{1}, 1), obj.maxIterations + 1);
            for i = 1:(obj.maxIterations + 1)
                Splt(:, i) = diag(obj.S{i});
            end
            p_S = plot(iterations, Splt, 'k');
            set(gca, 'fontsize', 15)
            set(p_S, 'linewidth', 1)
            axis square
            grid on
            xlabel('iteration i')
            ylabel('S')
            xlim([1 obj.maxIterations + 1])
            ax = gca;
            ax.XMinorGrid = 'on';
            ax.YMinorGrid = 'on';
            
            subplot(2, 2, 4)
            p_sw = semilogy(iterations, obj.MCMCStepWidth, 'k');
            set(gca, 'fontsize', 15)
            set(p_sw, 'linewidth', 1)
            axis square
            grid on
            xlabel('iteration i')
            ylabel('step width')
            xlim([1 obj.maxIterations + 1])
            ax = gca;
            ax.XMinorGrid = 'on';
            ax.YMinorGrid = 'on';
            
        end
        
    end
    
end
