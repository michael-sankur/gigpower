function [aPQ_noise, aZ_noise, aI_noise] = readNoise()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    noise = readmatrix('C:\Users\kathl\Desktop\LinDist3Flow\20180601\PYTHON\spu_noise.txt');

    % Read in noise matrices
    aPQ_noise = readmatrix('C:\Users\kathl\Desktop\LinDist3Flow\20180601\PYTHON\aPQ.txt');
    aPQ_noise = aPQ_noise(:, 1:end-1);

    aZ_noise = readmatrix('C:\Users\kathl\Desktop\LinDist3Flow\20180601\PYTHON\aZ.txt');
    aZ_noise = aZ_noise(:, 1:end-1);

    aI_noise = readmatrix('C:\Users\kathl\Desktop\LinDist3Flow\20180601\PYTHON\aI.txt');
    aI_noise = aI_noise(:, 1:end-1);
end

