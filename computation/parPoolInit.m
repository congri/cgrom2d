function [] = parPoolInit()
%% Initializes parallel pool

N_Threads = 16;
if isempty(gcp('nocreate'))
    % Create with N_Threads workers
    parpool('local', N_Threads);
end

end

