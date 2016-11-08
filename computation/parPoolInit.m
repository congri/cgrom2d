
function [] = parPoolInit(N_Threads)
%% Initializes parallel pool

if(nargin == 0 || N_Threads > 16)
    N_Threads = 16;
end
if isempty(gcp('nocreate'))
    % Create with N_Threads workers
    parpool('local', N_Threads);
end

end

