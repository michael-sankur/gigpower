function X = createfsolvefunction(feeder,nodes,lines,configs,loads,caps,controllers,type)

fid = fopen(['C:\Users\Michael\Desktop\reconfiguration\Mesh\TX\ftemp3phase' type '.m'],'w'); 

NPH = nodes.PH;
LPH = lines.PH;

nnode = nodes.nnode;
nline = lines.nline;

FZpu = lines.FZpu;

TXnum = lines.TXnum;
RXnum = lines.RXnum;

inmat = nodes.inmat;
outmat = nodes.outmat;

spu = loads.spu;
aPQ = loads.aPQ;
aZ = loads.aZ;

cappu = caps.cappu;

wpu = controllers.wpu;

%% Fsolve real

fprintf(fid,['function FT = ftemp3phase' type '(X)\n\n']);

fprintf(fid,['FT(1) = ' 'X(3) - 1' ';\n']);
fprintf(fid,['FT(2) = ' 'X(4) - 0' ';\n']);

fprintf(fid,['FT(3) = ' 'X(' num2str(2*nnode + 3) ') + 0.5' ';\n']);
fprintf(fid,['FT(4) = ' 'X(' num2str(2*nnode + 4) ') + sqrt(3)/2' ';\n']);

fprintf(fid,['FT(5) = ' 'X(' num2str(2*2*nnode + 3) ') + 0.5' ';\n']);
fprintf(fid,['FT(6) = ' 'X(' num2str(2*2*nnode + 4) ') - sqrt(3)/2' ';\n']);

%%

count = 7;
for ph = 1:3
    
    for k1 = 1:nline
        
%         idxre = 2*(ph-1)*nline + 2*k1-1;
%         idxim = 2*(ph-1)*nline + 2*k1;
        
        idxCmn = 2*3*nnode + 2*(ph-1)*nline + 2*k1-1;
        idxDmn = 2*3*nnode + 2*(ph-1)*nline + 2*k1;
        
        if LPH(ph,k1) == 0            
                        
            kvlre = ['X(' num2str(idxCmn) ');\n'];
%             fprintf(kvlre);
            fprintf(fid,['FT(' num2str(count) ') = ' kvlre]);
            count = count + 1;
            
            kvlim = ['X(' num2str(idxDmn) ');\n'];
%             fprintf(kvlim);
            fprintf(fid,['FT(' num2str(count) ') = ' kvlim]);
            count = count + 1;
        
        elseif LPH(ph,k1) == 1
                      
            idxAmTx = 2*(ph-1)*nnode + 2*TXnum(k1)-1;
            idxBmTx = 2*(ph-1)*nnode + 2*TXnum(k1);
            
            idxAmRx = 2*(ph-1)*nnode + 2*RXnum(k1)-1;            
            idxBmRx = 2*(ph-1)*nnode + 2*RXnum(k1);
            
            idxCmna = 2*3*nnode + 2*k1-1;
            idxDmna = 2*3*nnode + 2*k1;
            
            idxCmnb = 2*3*nnode + 2*nline + 2*k1-1;
            idxDmnb = 2*3*nnode + 2*nline + 2*k1;
            
            idxCmnc = 2*3*nnode + 2*2*nline + 2*k1-1;
            idxDmnc = 2*3*nnode + 2*2*nline + 2*k1;
          
            kvlre = ['X(' num2str(idxAmTx) ') - X(' num2str(idxAmRx) ')' ...
                ' - ' num2str(real(FZpu(ph,1,k1))) '*X(' num2str(idxCmna) ')' ...
                ' + ' num2str(imag(FZpu(ph,1,k1))) '*X(' num2str(idxDmna) ')' ...
                ' - ' num2str(real(FZpu(ph,2,k1))) '*X(' num2str(idxCmnb) ')' ...
                ' + ' num2str(imag(FZpu(ph,2,k1))) '*X(' num2str(idxDmnb) ')' ...
                ' - ' num2str(real(FZpu(ph,3,k1))) '*X(' num2str(idxCmnc) ')' ...
                ' + ' num2str(imag(FZpu(ph,3,k1))) '*X(' num2str(idxDmnc) ')' ...
                ';\n'];
        
%             fprintf(lineline);
            fprintf(fid,['FT(' num2str(count) ') = ' kvlre]);
            count = count + 1;
            
            kvlim = ['X(' num2str(idxBmTx) ') - X(' num2str(idxBmRx) ')' ...
                ' - ' num2str(real(FZpu(ph,1,k1))) '*X(' num2str(idxDmna) ')' ...
                ' - ' num2str(imag(FZpu(ph,1,k1))) '*X(' num2str(idxCmna) ')' ...
                ' - ' num2str(real(FZpu(ph,2,k1))) '*X(' num2str(idxDmnb) ')' ...
                ' - ' num2str(imag(FZpu(ph,2,k1))) '*X(' num2str(idxCmnb) ')' ...
                ' - ' num2str(real(FZpu(ph,3,k1))) '*X(' num2str(idxDmnc) ')' ...
                ' - ' num2str(imag(FZpu(ph,3,k1))) '*X(' num2str(idxCmnc) ')' ...
                ';\n'];
        
%             fprintf(lineline);
            fprintf(fid,['FT(' num2str(count) ') = ' kvlim]);
            count = count + 1;
                                
        end        
        
    end
        
end

%%

for ph = 1:3
    
    for k1 = 2:nnode
                
        idxAm = 2*(ph-1)*nnode + 2*k1-1;
        idxBm = 2*(ph-1)*nnode + 2*k1;
        
        if NPH(ph,k1) == 0
                      
            kclre = ['X(' num2str(idxAm) ');\n'];
%             fprintf([num2str(k1) ' - ' num2str(ph) ' - ' nodelineI]);
            fprintf(fid,['FT(' num2str(count) ') = ' kclre]);
            count = count + 1;
            
            kclim = ['X(' num2str(idxBm) ');\n'];
%             fprintf([num2str(k1) ' - ' num2str(ph) ' - ' nodelineI]);
            fprintf(fid,['FT(' num2str(count) ') = ' kclim]);
            count = count + 1;
            
        elseif NPH(ph,k1) == 1
            
            powinR = []; powinI = [];
            for k2 = 1:length(inmat(:,k1))
                
                idxClm = 2*3*nnode + 2*(ph-1)*nline + 2*inmat(k2,k1)-1;
                idxDlm = 2*3*nnode + 2*(ph-1)*nline + 2*inmat(k2,k1);
                
                if inmat(k2,k1) ~= 0 && k2 == 1                                        
                    powinR = ['X(' num2str(idxAm) ')*X(' num2str(idxClm) ') + X(' num2str(idxBm) ')*X(' num2str(idxDlm) ')'];                    
                    powinI = ['-X(' num2str(idxAm) ')*X(' num2str(idxDlm) ') + X(' num2str(idxBm) ')*X(' num2str(idxClm) ')'];
                end
                if inmat(k2,k1) ~= 0 && k2 >= 2                                        
                    powinR = [powinR ' + X(' num2str(idxAm) ')*X(' num2str(idxClm) ') + X(' num2str(idxBm) ')*X(' num2str(idxDlm) ')'];                    
                    powinI = [powinI '- X(' num2str(idxAm) ')*X(' num2str(idxDlm) ') + X(' num2str(idxBm) ')*X(' num2str(idxClm) ')'];
                end
                
            end
%             fprintf([num2str(k1) ' - ' num2str(ph) ' - ' powinR '\n'])
%             fprintf([num2str(k1) ' - ' num2str(ph) ' - ' powinI '\n'])
            
            sR = [' - ' num2str(real(spu(ph,k1))) '*(' num2str(aPQ(ph,k1)) ' + ' ...
                num2str(aZ(ph,k1)) '*(X(' num2str(idxAm) ')^2 + X(' num2str(idxBm) ')^2))' ...
                ' - ' num2str(real(wpu(ph,k1)))];
            
            sI = [' - ' num2str(imag(spu(ph,k1))) '*(' num2str(aPQ(ph,k1)) ' + ' ...
                num2str(aZ(ph,k1)) '*(X(' num2str(idxAm) ')^2 + X(' num2str(idxBm) ')^2))' ...
                ' - ' num2str(imag(wpu(ph,k1))) ' + ' num2str(cappu(ph,k1))];
                
                
%             fprintf([num2str(k1) ' - ' num2str(ph) ' - ' sR '\n'])
%             fprintf([num2str(k1) ' - ' num2str(ph) ' - ' sI '\n'])
            
            powoutR = []; powoutI = [];
            for k2 = 1:length(outmat(:,k1))
                
                idxCmn = 2*3*nnode + 2*(ph-1)*nline + 2*outmat(k2,k1)-1;
                idxDmn = 2*3*nnode + 2*(ph-1)*nline + 2*outmat(k2,k1);
                
                if outmat(k2,k1) ~= 0 && k2 == 1
                    powoutR = [' - X(' num2str(idxAm) ')*X(' num2str(idxCmn) ') - X(' num2str(idxBm) ')*X(' num2str(idxDmn) ')'];                    
                    powoutI = [' + X(' num2str(idxAm) ')*X(' num2str(idxDmn) ') - X(' num2str(idxBm) ')*X(' num2str(idxCmn) ')'];
                end
                if outmat(k2,k1) ~= 0 && k2 >= 2
                    powoutR = [powoutR ' - X(' num2str(idxAm) ')*X(' num2str(idxCmn) ') - X(' num2str(idxBm) ')*X(' num2str(idxDmn) ')'];                    
                    powoutI = [powoutI ' + X(' num2str(idxAm) ')*X(' num2str(idxDmn) ') - X(' num2str(idxBm) ')*X(' num2str(idxCmn) ')'];
                end
                
            end
%             fprintf([num2str(k1) ' - ' num2str(ph) ' - ' powoutR '\n'])
%             fprintf([num2str(k1) ' - ' num2str(ph) ' - ' powoutI '\n'])
            
            kclre = [powinR sR powoutR ';\n'];
%             fprintf([num2str(k1) ' - ' num2str(ph) ' - ' kclre]);
            fprintf(fid,['FT(' num2str(count) ') = ' kclre]);
            count = count + 1;
            
            
            kclim = [powinI sI powoutI ';\n'];
%             fprintf([num2str(k1) ' - ' num2str(ph) ' - ' kclim]);
            fprintf(fid,['FT(' num2str(count) ') = ' kclim]);
            count = count + 1;
            
        end
    end
end

fclose(fid);

end