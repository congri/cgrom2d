classdef FinescaleData
    %This class holds all the properties of the finescale data set and the methods to generate
    %the finescale data
    
    properties
        distributionType = 'correlated_binary';     %Type of p(x); usually use 'correlated_binary', which samples
                                                    %the conductivity field from a spacially
                                                    %correlated Gaussian
        nSets                                       %Number of data sets to generate;
                                                    %only difference is that different sets are
                                                    %stored to different .mat files
        nSamples                                    %Vector giving the number of samples per set
    end
    
    methods
    end
    
end

