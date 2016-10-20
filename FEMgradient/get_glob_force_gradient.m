function [F] = get_glob_force_gradient(domain, k)
%Assemble global force vector

%nested function for performance
f2 = zeros(4, domain.nEl);
Tb = zeros(4, 1);
    function f2 = get_loc_force_gradient2(domain, k)
        f2 = zeros(4, domain.nEl);
        for ee = 1:domain.nEl
            %Boundary value temperature of element e
            Tb = zeros(4, 1);
            Tbflag = false;
            for i = 1:4
                globNode = domain.globalNodeNumber(ee, i);
                if(~isnan(domain.essentialTemperatures(globNode)))
                    Tb(i) = domain.essentialTemperatures(globNode);
                    Tbflag = true;
                end
            end

            if(Tbflag)
                f2(:, ee) = -k(:, :, ee)*Tb;
            end
        end
    end


F = zeros(domain.nEq,1);

f = get_loc_force_gradient2(domain, k);
for e = 1:domain.nEl
    for ln = 1:4
       eqn = domain.lm(e,ln);
       if(eqn)
          F(eqn) = F(eqn) + f(ln,e);
       end
    end
end

end