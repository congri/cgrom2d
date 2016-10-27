%Generate boundary condition functions

%% Temperature field and gradient generating the boundary conditions
a = [-5 3 -2 0];
Tb = @(x) a(1) + a(2)*x(1) + a(3)*x(2) + a(4)*x(1)*x(2);
qb{1} = @(x) -(a(3) + a(4)*x);      %lower bound
qb{2} = @(y) (a(2) + a(4)*y);       %right bound
qb{3} = @(x) (a(3) + a(4)*x);       %upper bound
qb{4} = @(y) -(a(2) + a(4)*y);      %left bound