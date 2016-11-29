function [curr_x, grad, Hess, nIter] = newtonRaphsonMaximization(objective, startValue, Xtol, provide_objective, debug)
%Maximization based on the well-known Newton-Raphson method
%%Input:
%   objective:              func. handle with [grad., Hess., obj.] as output
%   startValue:             Optimization start value
%   Xtol:                   Tolerance in x for convergence

%step size parameter, 0 < gamma < 1
gamma = .5;
RMfac2 = (1e-5)*ones(length(startValue), 1);

%%Optimization
converged = false;
curr_x = startValue;        %clear notation
if provide_objective
    [grad, Hess, obj] = objective(curr_x);
else
    [grad, Hess] = objective(curr_x);
end
if(debug)
    gradArray = grad;
end
step = Hess\grad;
nIter = 1;
%Iterate until convergence
while(~converged)
    old_x = curr_x;
    if(~provide_objective)
        curr_x = old_x - gamma*step;
    else
        temp_x = old_x - gamma*step;
        [grad_temp, Hess_temp, obj_temp] = objective(temp_x);
        nIter = nIter + 1;
        step_temp = Hess_temp\grad_temp;
        while(obj_temp < obj)
            step_temp = .9*step_temp;
            temp_x = old_x - gamma*step_temp;
            [grad_temp, Hess_temp, obj_temp] = objective(temp_x);
            nIter = nIter + 1;
        end
        %step accepted
        curr_x = temp_x;
        obj = obj_temp;
        grad = grad_temp;
        Hess = Hess_temp;
        step = Hess\grad;
    end
    
    %check for convergence
    if(norm(curr_x - old_x)/norm(curr_x) < Xtol)
        converged = true;
    else    %if not converged, compute next step
        if(~provide_objective)
            grad_old = grad;
            [grad, Hess] = objective(curr_x);
            %check for negative definiteness
            if(~all(all(isfinite(Hess))) || any(eig(Hess) >= 0))
                warning('Hessian not negative definite, do Robbins-Monro instead')
                RMoff = 100;
                RMfac1 = RMoff/(RMoff + nIter);
                signChange = grad_old.*grad;
                signChange = (signChange < 0);
                for j = 1:length(grad)
                   if(signChange(j))
                       %decrease stepsize when gradient sign changes
                       RMfac2(j) = .9*RMfac2(j); 
                   else
                       %increase step size when gradient sign does not change
                       RMfac2(j) = 1.1*RMfac2(j);
                   end
                end
                RMstep = RMfac1*RMfac2;
                step = - norm(curr_x)*(RMstep.*(grad/norm(grad)));
            else
                step = Hess\grad;
            end
            nIter = nIter + 1;
        end
        if(debug)
            curr_x
            grad
            Hess
            invHess = inv(Hess)
            step
            gradArray = [gradArray, grad];
            plot(gradArray')
            xlim([1 nIter + 5])
            drawnow
            if(mod(nIter, 500) == 0)
                pause
            end
        end
    end
end

end

