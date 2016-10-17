function [cond, Tf] = genData(domain, physical, fineData)
%Generating full order data

%% Draw conductivity/ log conductivity
if strcmp(fineData.dist, 'uniform')
    %conductivity uniformly distributed between lo and up
    cond = (fineData.up - fineData.lo)*rand(domain.nEl, fineData.nSamples) + fineData.lo;
elseif strcmp(fineData.dist, 'gaussian')
    %log conductivity gaussian distributed
    x = normrnd(fineData.mu, fineData.sigma, domain.nEl, fineData.nSamples);
    cond = exp(x);
elseif strcmp(fineData.dist, 'binary')
    %binary distribution of conductivity (Bernoulli)
    r = rand(domain.nEl, fineData.nSamples);
    cond = fineData.lo*ones(domain.nEl, fineData.nSamples);
    cond(r > fineData.p_lo) = fineData.up;
elseif strcmp(fineData.dist, 'predefined_binary')
    cond(:, 1) = fineData.lo*ones(domain.nEl, 1);
    cond(:, 2) = fineData.up*ones(domain.nEl, 1);
    %Define blocks. Only works for nf = 16!
    block1 = [18 19 20 34 35 36 50 51 52];
    block2 = block1 + 7;
    block3 = [186 187 188 189 201 202 203 215 216 217 232];
    block4 = [150 151 152 153 154 155 166 167 168 169 170];
    block5 = [18 19 20 33 34 35 36 37 38 39 40 50 51 52 67 83];
    %cross in lower left badges
    block6 = [33:40, 4 20 36 52 68 84 100 116];
    %checkerboard in lower right
    block7 = [9 11 13 15, 26 28 30 32, 41 43 45 47, 58 60 62 64,...
                73 75 77 79, 90 92 94 96, 105 107 109 111, 122 124 126 128];
    %diagonal cross in upper left
    block8 = [241   224   207   190   173   156   139   122, 129:17:248];
    %center square in upper right
    block9 = [154:159 170:175 186:191 202:207 218:223 234:239];
    
    %Matrix phases
    Matrix1 = fineData.lo*ones(domain.nEl, 1);
    Matrix2 = fineData.up*ones(domain.nEl, 1);
    %All blocks in matrix phase
    cond(:, 3) = Matrix1;
    cond([block1, block2, block3, block4], 3) = fineData.up;
    cond(:, 4) = Matrix2;
    cond([block1, block2, block3, block4], 4) = fineData.lo;
    cond(:, 5) = Matrix1;
    cond([block1, block4], 5) = fineData.up;
    cond(:, 6) = Matrix2;
    cond([block1, block4], 6) = fineData.lo;
    cond(:, 7) = Matrix1;
    cond([block1 block2 block3 block4 block5], 7) = fineData.up;
    cond(:, 8) = Matrix2;
    cond([block1 block2 block3 block4 block5], 8) = fineData.lo;
    cond(:, 9) = Matrix1;
    cond([block6 block7 block8 block9], 9) = fineData.up;
    cond(:, 10) = Matrix2;
    cond([block6 block7 block8 block9], 10) = fineData.lo;
    
    
    %The rest is random
    r = rand(domain.nEl, 6);
    condtemp = fineData.lo*ones(domain.nEl, 6);
    condtemp(r > .5) = fineData.up;
    cond(:, 11:16) = condtemp;
else
    error('unknown FOM conductivity distribution');
end

%% Compute output data (finescale nodal temperatures)
Tf = zeros(domain.nNodes, fineData.nSamples);
D = zeros(2, 2, domain.nEl);
for i = 1:fineData.nSamples
    control.plt = false;
    %Conductivity matrix D, only consider isotropic materials here
    for j = 1:domain.nEl
        D(:, :, j) =  cond(j, i)*eye(2);
    end
    FEMout = heat2d(domain, physical, control, D);
    %Store fine temperatures as a vector Tf. Use reshape(Tf(:, i), domain.nElX + 1, domain.nElY + 1)
    %and then transpose result to reconvert it to original temperature field
    Ttemp = FEMout.Tff';
    Tf(:, i) = Ttemp(:);
end

end

