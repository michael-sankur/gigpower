function [VNR, INR, STXNR, SRXNR, iNR, sNR, vvcNR, vvciter] = NR3VVC(network1,Vopt,Iopt,slackidx,slackVnom)

nnode = network1.nodes.nnode;

VNRVVC0 = ones(3,nnode);
VNR = zeros(3,nnode);

vvciter = 0;
while max(max(abs(VNRVVC0 - VNR))) >= 1e-12
        
    VNRVVC0 = VNR;
    
    network1.vvc.vvcpu = zeros(3,nnode);
    for ph = 1:3
        for kn = 2:nnode
            if network1.vvc.state(ph,kn) == 1
                qk = VVC_corrected(abs(VNRVVC0(ph,kn)),network1.vvc.qminpu(ph,kn),network1.vvc.qmaxpu(ph,kn),network1.vvc.Vmin(ph,kn),network1.vvc.Vmax(ph,kn));
                network1.vvc.vvcpu(ph,kn) = qk;
            end        
        end
    end    
    network1.vvc.vvcpu;
    
    [VNR, INR, STXNR, SRXNR, iNR, sNR, nriter] = NR3(network1,Vopt,Iopt,slackidx,slackVnom);
    vvciter = vvciter + 1;
    
end

vvcNR = network1.vvc.vvcpu;

end