% Michael Sankur - msankur@lbl.gov
% 2016.06.03

% This is a one time step simulation in which DERs are dispatched to
% balance the voltage across the existing phases on each node of the
% feeder.

clc, clear all, close all

if exist('FBS','file') == 0 || exist('feeder_mapper_ieee_function_20160603','file') == 0
    path(path,'C:\Users\Michael\Desktop\reconfiguration\20160603');
end

if exist('mxlpsolve','file') == 0
    path(path,'C:\Users\Michael\Desktop\mxlp');
end

if exist('cvx_begin','file') == 0
    cd C:\Users\Michael\Desktop\cvx-w64\cvx
    cvx_setup
end

%% Load feeder

% name = '6node_fullphase_test';
name = 'ieee_13node_new';

fp = ['C:\Users\Michael\Desktop\reconfiguration\Feeders\'];
fn = [name '.txt'];

[feeder, nodes, lines, configs, loads, caps, controllers] = feeder_mapper_TX_function_20170215(name, fp, fn);


%% Feeder paramaters

nnode = nodes.nnode;
nline = lines.nline;

%% Load parameters

loads.aPQ = 0.85*ones(3,nnode).*nodes.PH;
loads.aI = zeros(3,nnode);
loads.aZ = 0.15*ones(3,nnode).*nodes.PH;

loads.spu = 1.0*loads.spu;

%% Capacitor parameters

caps.cappu = 0*caps.cappu;

%% Controller parameters

controllers.wmaxpu = 1*controllers.wmaxpu;

controllers.wpu = zeros(3,nnode);

%% FBS New

Vnom = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];

%% Monte Carlo

Plim = 0.00:0.01:0.15;
Qlim = 0.00:0.01:0.15;

kp = 7;
kq = 8;

SS = [];
ME = [];
TE = [];
VE = [];
PE = [];
QE = [];
SE = [];

tic

for kp = 2:length(Plim)
    for kq = 2:length(Qlim)
        
%         kp = length(Plim);
%         kq = length(Qlim);
        
        for kmc = 1:1

            for knode = 3:nnode
                for ph = 1:3

                    if nodes.PH(ph,knode) == 1
                        loads.spu(ph,knode) = rand*Plim(kp) + 1j*rand*Qlim(kq);
                    end
                    
                    if rand < 0.1
                        loads.spu(ph,knode) = 0;
                    end

                end
            end

            FBS0 = FBS_RX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
            FBS1 = FBS_RX_A1(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
            
            SS(end+1) = sum(abs(FBS0.Srx(:,2)));
            
            [a, b] = max(abs(abs(FBS0.V(:)) - abs(FBS1.V(:))));
            ME(end+1) = a/abs(FBS0.V(b));
            
            TE(end+1) = max(max(abs(angle(FBS0.V) - angle(FBS1.V))));
            
            [a, b] = max(abs(FBS0.V(:) - FBS1.V(:)));
            VE(end+1) = a/abs(FBS0.V(b));
            
            [a, b] = max(abs(real(FBS0.Srx(:)) - real(FBS1.Srx(:))));
            PE(end+1) = a/abs(real(FBS0.Srx(b)));
            
            [a, b] = max(abs(imag(FBS0.Srx(:)) - imag(FBS1.Srx(:))));
            QE(end+1) = a/abs(imag(FBS1.Srx(b)));
            
            [a, b] = max(abs(abs(FBS0.Srx(:)) - abs(FBS1.Srx(:))));
            SE(end+1) = a/abs(FBS0.Srx(b));
            
        end
        
    end
end

toc

%%

close all

figure,
plot(SS,ME,'r.','LineWidth',2)
set(gca,'YTick',0:0.005:0.02,'Fontsize',15,'FontWeight','bold')
title('Maxiumum Voltage Magnitude Error','Fontsize',15,'FontWeight','bold')
xlabel('Substation power |S_{0}|','Fontsize',15,'FontWeight','bold')
ylabel('Magnitude error [p.u.]','Fontsize',15,'FontWeight','bold')
axis([0 1.5 0 0.02])
pbaspect([1 0.25 1])

print('-f1','-depsc',['C:\Users\Michael\Desktop\temp\montecarlomagnitude.eps'])


figure,
plot(SS,180/pi*TE,'r.','LineWidth',2)
set(gca,'YTick',0:0.25:1,'Fontsize',15,'FontWeight','bold')
title('Maxiumum Voltage Angle Error','Fontsize',15,'FontWeight','bold')
xlabel('Substation power |S_{0}|','Fontsize',15,'FontWeight','bold')
ylabel('Angle error [deg]','Fontsize',15,'FontWeight','bold')
axis([0 1.5 0 1])
pbaspect([1 0.25 1])

print('-f2','-depsc',['C:\Users\Michael\Desktop\temp\montecarloangle.eps'])


figure,
plot(SS,VE,'r.','LineWidth',2)

% figure,
% plot(SS,PE,'r.','LineWidth',2)
% 
% figure,
% plot(SS,QE,'r.','LineWidth',2)

figure,
plot(SS,SE,'r.','LineWidth',2)
set(gca,'YTick',0:0.025:0.1,'Fontsize',15,'FontWeight','bold')
title('Maxiumum Apparent Power Error','Fontsize',15,'FontWeight','bold')
xlabel('Substation power |S_{0}|','Fontsize',15,'FontWeight','bold')
ylabel('Power error [p.u.]','Fontsize',15,'FontWeight','bold')
axis([0 1.5 0 0.1])
pbaspect([1 0.25 1])

print('-f3','-depsc',['C:\Users\Michael\Desktop\temp\montecarlopower.eps'])



% lines.origFZpu = lines.FZpu;
% loads.origspu = loads.spu;
% 
% ma1 = 0.5:0.025:2;
% ma2 = 0.75:0.025:1.25;
% 
% EERROR = [];
% DERROR = [];
% VERROR = [];
% PERROR = [];
% QERROR = [];
% SERROR = [];
% 
% for k1 = 1:length(ma1)
%     for k2 = 1:length(ma2)
%     
%         loads.spu = ma1(k1)*loads.origspu;
%         lines.FZpu = ma2(k2)*lines.origFZpu;
% 
%         FBS0 = FBS_RX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% 
%         FBS1 = FBS_RX_A1(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% 
% %         FBS2 = FBS_RX_A2(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% % 
% %         FBS3 = FBS_RX_A3(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% % 
% %         FBS4 = FBS_RX_A4(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% % 
% %         FBS5 = FBS_RX_A5(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% % 
% %         FBS6 = FBS_RX_A6(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
% 
%         EERROR(k1,k2) = max(max(abs(abs(FBS0.V) - abs(FBS1.V))));
%         DERROR(k1,k2) = max(max(abs(angle(FBS0.V) - angle(FBS1.V))));
%         VERROR(k1,k2) = max(max(abs(FBS0.V - FBS1.V)));
%         
%         PERROR(k1,k2) = max(max(abs(real(FBS0.Srx - FBS1.Srx))));
%         QERROR(k1,k2) = max(max(abs(imag(FBS0.Srx - FBS1.Srx))));
%         SERROR(k1,k2) = max(max(abs(FBS0.Srx - FBS1.Srx)));
% 
%     end
% end
% 
% %%
% 
% lss = {'r.','g.','b.','m.','k.','c.'};
% 
% figure, box on, hold on
% % surf(ma1,ma2,EERROR','LineStyle','none')
% surf(ma1,ma2,DERROR','LineStyle','none')
