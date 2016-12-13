function [PhiArray] = normalizeDesignMatrix(PhiArray, colNormPhi)
%Divides each column in PhiArray by the value stored in colNormPhi


for i = 1:size(PhiArray, 2)
    for j = 1:size(PhiArray, 3)
        PhiArray(:, i, j) = PhiArray(:, i, j)/colNormPhi(:, i);
    end
end


end

