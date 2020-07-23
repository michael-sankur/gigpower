function [network1] = addNoise(network1, kt, nnode, aPQ, aZ, aI)
%Adds noise to aPQ, aZ, aI loads in network 

    for i = 1:3
        for j = 1:nnode
            network1.loads.aPQ = aPQ((3*kt-2):kt*3, :) .*(network1.loads.spu ~= 0);
            network1.loads.aZ = aZ((3*kt-2):kt*3, :) .*(network1.loads.spu ~= 0);
            network1.loads.aI = aI((3*kt-2):kt*3, :) .*(network1.loads.spu ~= 0);
        end
    end
   
end

