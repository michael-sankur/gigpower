function [aPQ, aZ, aI, spu] = readNoise()
%Reads spu, aPQ, aZ, aI noise in from text files generated in
%NR3test.ipynb


    spu = readmatrix('C:\Users\kathl\Desktop\LinDist3Flow\20180601\PYTHON\spu_noise.txt');

    % Read in noise matrices
    aPQ = readmatrix('C:\Users\kathl\Desktop\LinDist3Flow\20180601\PYTHON\aPQ.txt');
    aPQ = aPQ(:, 1:end-1);

    aZ = readmatrix('C:\Users\kathl\Desktop\LinDist3Flow\20180601\PYTHON\aZ.txt');
    aZ = aZ(:, 1:end-1);

    aI = readmatrix('C:\Users\kathl\Desktop\LinDist3Flow\20180601\PYTHON\aI.txt');
    aI = aI(:, 1:end-1);
end

