function [Nbubbles] = numberOfObjects(lambdaMat, fineData, hilo)
%Counts the number of disconnected conductivity bubbles
%   lambdaMat:     2D(!) conductivity image, i.e. 2-dim array

if strcmp(hilo, 'hi')
    cc = bwconncomp((lambdaMat > fineData.lo));
elseif strcmp(hilo, 'lo')
    cc = bwconncomp((lambdaMat < fineData.up));
else
    error('High or low conductivity phase bubbles?')
end

Nbubbles = cc.NumObjects;

end

