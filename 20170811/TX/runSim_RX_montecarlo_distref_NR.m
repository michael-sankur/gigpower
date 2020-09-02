% Michael Sankur - msankur@lbl.gov
% 2016.06.03

% This is a one time step simulation in which DERs are dispatched to
% balance the voltage across the existing phases on each node of the
% feeder.

clc, clear all, close all

% if exist('FBS','file') == 0 || exist('feeder_mapper_ieee_function_20160603','file') == 0
%     path(path,'C:\Users\Michael\Desktop\reconfiguration\20160603');
% end

path(path,genpath('C:\Users\Michael\Desktop\reconfiguration\Mesh\GEN'));
path(path,'C:\Users\Michael\Desktop\reconfiguration\Mesh\SDP');
path(path,'C:\Users\Michael\Desktop\reconfiguration\Mesh\TX');
% path(path,'C:\Users\Michael\Desktop\heatscatter');

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
% name = 'ieee_37node_new';
% name = 'mesh_09node_balance';

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

caps.cappu = 1*caps.cappu;

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
% MEM = []; ME2 = [];
% EEM = []; EE2 = [];
% DEM = []; DE2 = [];
% TEM = []; TE2 = [];
% PEM = []; PE2 = [];
% QEM = []; QE2 = [];
% SEM = []; SE2 = [];
% 
% EEM2 = []; EE22 = [];
% DEM2 = []; DE22 = [];
% SEM2 = []; SE22 = [];

EEM1 = []; EEM1P = [];
DEM1 = []; DEM1P = [];
DEM2 = []; DEM2P = [];

PEM1 = []; PEM1P = [];
QEM1 = []; QEM1P = [];
SEM1 = []; SEM1P = [];

TVEM1 = [];
TVEM2 = [];

VV = [];

tic

for kp = 2:length(Plim)
%     kp
    for kq = 2:length(Qlim)
%         kq
        
%         kp = length(Plim);
%         kq = length(Qlim);
        
        for kmc = 1:25
%             kmc

            for knode = 3:nnode
                for ph = 1:3

                    if nodes.PH(ph,knode) == 1
%                         loads.spu(ph,knode) = rand*Plim(kp) + 1j*rand*Qlim(kq);
                        loads.spu(ph,knode) = (0.5*Plim(kp) + randn*0.25*Plim(kp)) + 1j*(0.5*Qlim(kq) + 0.25*randn*Qlim(kq));
                    end
                    
                    if rand < 0.1
                        loads.spu(ph,knode) = 0;
                    end

                end
            end

%             FBS0 = FBS_RX(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
%             FBS1 = FBS_RX_A1(feeder,nodes,lines,configs,loads,caps,controllers,Vnom);
            
            [VNR0, INR0, STXNR0, SRXNR0] = NR3(feeder,nodes,lines,configs,loads,caps,controllers,[],[],1,Vnom);
            
            [VNRA1, STXNRA1, SRXNRA1] = NR3A1(feeder,nodes,lines,configs,loads,caps,controllers,[],[],1,Vnom);
            
            emeas = abs(VNR0);
            [VNRA2, STXNRA2, SRXNRA2] = NR3A2(feeder,nodes,lines,configs,loads,caps,controllers,[],[],emeas,1,Vnom);
            
            SS(end+1) = sum(abs(STXNR0(:,1)));
            
%             [a, b] = max(abs(abs(VNR0(:)) - abs(VNRA1(:))));
%             MEM(end+1) = a/abs(VNR0(b));
%             
%             TEM(end+1) = max(max(abs(angle(VNR0) - angle(VNRA1))));
%             
%             [a, b] = max(abs(VNR0(:) - VNRA1(:)));
%             VEM(end+1) = a/abs(VNR0(b));
%             
%             [a, b] = max(abs(real(SRXNR0(:)) - real(SRXNRA1(:))));
%             PEM(end+1) = a/abs(real(SRXNR0(b)));
%             
%             [a, b] = max(abs(imag(SRXNR0(:)) - imag(SRXNRA1(:))));
%             QEM(end+1) = a/abs(imag(SRXNR0(b)));
%             
%             [a, b] = max(abs(abs(SRXNR0(:)) - abs(SRXNRA1(:))));
%             SEM(end+1) = a/abs(SRXNR0(b));
            
            VNR0 = [VNR0(1,:);
                exp(j*120*pi/180)*VNR0(2,:);
                exp(j*240*pi/180)*VNR0(3,:)];
            
            VNRA1 = [VNRA1(1,:);
                exp(j*120*pi/180)*VNRA1(2,:);
                exp(j*240*pi/180)*VNRA1(3,:)];
            
            VNRA2 = [VNRA2(1,:);
                exp(j*120*pi/180)*VNRA2(2,:);
                exp(j*240*pi/180)*VNRA2(3,:)];
            
            if sum(abs(VNR0(:)) <= 0.94) > sum(abs(VNR0(:)) == 0)
                VV(end+1) = 1;
                abs(VNR0);
            elseif sum(abs(VNR0(:)) <= 0.94) == sum(abs(VNR0(:)) == 0)
                VV(end+1) = 0;
            end
            
            
            [a, b] = max(abs(abs(VNR0(:)) - abs(VNRA1(:))));
            
            EEM1(end+1) = a;
%             EEM1P(end+1) = abs(a/(abs(VNR0(b)) - 1));
            
            [a, b] = max(abs(angle(VNR0(:)) - angle(VNRA1(:))));
            
            DEM1(end+1) = a;
%             DEM1P(end+1) = abs(a/angle(VNR0(b)));
            
            [a, b] = max(abs(angle(VNR0(:)) - angle(VNRA2(:))));
            
            DEM2(end+1) = a;
%             DEM2P(end+1) = abs(a/angle(VNR0(b)));

            [a, b] = max(abs(VNR0(:) - VNRA1(:)));

            TVEM1(end+1) = a;
            
            [a, b] = max(abs(VNR0(:) - VNRA2(:)));

            TVEM2(end+1) = a;
            
            
            
                        
            [a, b] = max(abs(real(SRXNR0(:)) - real(SRXNRA1(:))));
            
            PEM1(end+1) = a;
%             PEM1P(end+1) = abs(a/real(SRXNR0(b)));
            
            [a, b] = max(abs(imag(SRXNR0(:)) - imag(SRXNRA1(:))));
            
            QEM1(end+1) = a;
%             QEM1P(end+1) = abs(a/imag(SRXNR0(b)));
            
            [a, b] = max(abs(SRXNR0(:) - SRXNRA1(:)));
            
            SEM1(end+1) = a;
%             SEM1P(end+1) = abs(a/abs(SRXNR0(b)));
            
            
            
%             EEM1(end+1) = norm(abs(VNR0(:)) - abs(VNRA1(:)),Inf);
%             EEM2(end+1) = norm(abs(VNR0(:)) - abs(VNRA2(:)),Inf);
%             
%             
%             DEM1(end+1) = norm(angle(VNR0(:)) - angle(VNRA1(:)),Inf);            
%             DEM2(end+1) = norm(angle(VNR0(:)) - angle(VNRA2(:)),Inf);
%             
%             
%             PEM1(end+1) = norm(real(SRXNR0(:)) - real(SRXNRA1(:)),Inf);            
%             PEM2(end+1) = norm(real(SRXNR0(:)) - real(SRXNRA1(:)),Inf);
%             
%             
%             QEM1(end+1) = norm(imag(SRXNR0(:)) - imag(SRXNRA1(:)),Inf);            
%             QEM2(end+1) = norm(imag(SRXNR0(:)) - imag(SRXNRA1(:)),Inf);
%             
%             
%             SEM1(end+1) = norm(abs(SRXNR0(:)) - abs(SRXNRA1(:)),Inf);            
%             SEM2(end+1) = norm(abs(SRXNR0(:)) - abs(SRXNRA1(:)),Inf);
            

%             EEM(end+1) = norm(abs(VNR0(:)) - abs(VNRA1(:)),Inf);
%             EE2(end+1) = norm(abs(VNR0(:)) - abs(VNRA1(:)),2);
%             
%             DEM(end+1) = norm(angle(VNR0(:)) - angle(VNRA1(:)),Inf);
%             DE2(end+1) = norm(angle(VNR0(:)) - angle(VNRA1(:)),2);
%             
%             TEM(end+1) = norm(VNR0(:) - VNRA1(:),Inf);
%             TE2(end+1) = norm(VNR0(:) - VNRA1(:),2);
%             
%             PEM(end+1) = norm(real(SRXNR0(:)) - real(SRXNRA1(:)),Inf);
%             PE2(end+1) = norm(real(SRXNR0(:)) - real(SRXNRA1(:)),2);
%             
%             QEM(end+1) = norm(imag(SRXNR0(:)) - imag(SRXNRA1(:)),Inf);
%             QE2(end+1) = norm(imag(SRXNR0(:)) - imag(SRXNRA1(:)),2);
%             
%             SEM(end+1) = norm(abs(SRXNR0(:)) - abs(SRXNRA1(:)),Inf);
%             SE2(end+1) = norm(abs(SRXNR0(:)) - abs(SRXNRA1(:)),2);
%             
%             
%             EEM2(end+1) = norm(abs(VNR0(:)) - abs(VNRA2(:)),Inf);
%             EE22(end+1) = norm(abs(VNR0(:)) - abs(VNRA2(:)),2);
%             
%             DEM2(end+1) = norm(angle(VNR0(:)) - angle(VNRA2(:)),Inf);
%             DE22(end+1) = norm(angle(VNR0(:)) - angle(VNRA2(:)),2);
% %             
% %             TEM(end+1) = norm(VNR0(:) - VNRA1(:),Inf);
% %             TE2(end+1) = norm(VNR0(:) - VNRA1(:),2);
% %             
% %             PEM(end+1) = norm(real(SRXNR0(:)) - real(SRXNRA1(:)),Inf);
% %             PE2(end+1) = norm(real(SRXNR0(:)) - real(SRXNRA1(:)),2);
% %             
% %             QEM(end+1) = norm(imag(SRXNR0(:)) - imag(SRXNRA1(:)),Inf);
% %             QE2(end+1) = norm(imag(SRXNR0(:)) - imag(SRXNRA1(:)),2);
% %             
%             SEM2(end+1) = norm(abs(SRXNR0(:)) - abs(SRXNRA2(:)),Inf);
%             SE22(end+1) = norm(abs(SRXNR0(:)) - abs(SRXNRA2(:)),2);
            
            
        end
        
    end
end

toc

save('C:\Users\Michael\Desktop\temp\MCRES.mat','SS','EEM1','DEM1','SEM1')


%%

close all

if exist('EEM1','var') == 0 && exist('DEM1','var') == 0 && exist('SEM1','var') == 0
    load('C:\Users\Michael\Desktop\temp\MCRES.mat')
end

% figure,
% % plot(SS(VV == 0),EEM1(VV == 0),'g.',SS(VV == 1),EEM1(VV == 1),'r.','LineWidth',2)
% plot(SS,EEM1,'r.','LineWidth',2)
% set(gca,'YTick',0:0.005:0.02,'Fontsize',15,'FontWeight','bold')
% title('Maxiumum Voltage Magnitude Error','Fontsize',15,'FontWeight','bold')
% xlabel('Substation power |S_{0}|','Fontsize',15,'FontWeight','bold')
% % ylabel('Magnitude error [p.u.]','Fontsize',15,'FontWeight','bold')
% ylabel([char(949) '_{mag} [p.u.]'],'Fontsize',15,'FontWeight','bold')
% axis([0 1.5 0 0.015])
% pbaspect([1 0.35 1])
% 
% print('-f1','-depsc',['C:\Users\Michael\Desktop\temp\mcmagerror.eps'])
% print('-f1','-dpng',['C:\Users\Michael\Desktop\temp\mcmagerror.png'])
% 
% 
% % figure,
% % plot(SS(VV == 0),EEM1P(VV == 0),'g.',SS(VV == 1),EEM1P(VV == 1),'r.','LineWidth',2)
% % set(gca,'Fontsize',15,'FontWeight','bold')
% % title('Maxiumum Voltage Magnitude Error','Fontsize',15,'FontWeight','bold')
% % xlabel('Substation power |S_{0}|','Fontsize',15,'FontWeight','bold')
% % ylabel('Magnitude error [p.u.]','Fontsize',15,'FontWeight','bold')
% % axis([0 1.5 0 inf])
% % pbaspect([1 0.25 1])
% 
% 
% figure,
% % plot(SS(VV == 0),180/pi*DEM1(VV == 0),'g.',SS(VV == 1),180/pi*DEM1(VV == 1),'r.','LineWidth',2)
% plot(SS,180/pi*DEM1,'g.','LineWidth',2)
% set(gca,'YTick',0:0.1:0.5,'Fontsize',15,'FontWeight','bold')
% title('Maxiumum Voltage Angle Error','Fontsize',15,'FontWeight','bold')
% xlabel('Substation power |S_{0}|','Fontsize',15,'FontWeight','bold')
% % ylabel('Angle error [deg]','Fontsize',15,'FontWeight','bold')
% ylabel([char(949) '_{angle} [p.u.]'],'Fontsize',15,'FontWeight','bold')
% axis([0 1.5 0 0.5])
% pbaspect([1 0.375 1])
% 
% print('-f2','-depsc',['C:\Users\Michael\Desktop\temp\mcangleerror1.eps'])
% print('-f2','-dpng',['C:\Users\Michael\Desktop\temp\mcangleerror1.png'])
% 
% % figure,
% % plot(SS(VV == 0),DEM1P(VV == 0),'g.',SS(VV == 1),DEM1P(VV == 1),'r.','LineWidth',2)
% % set(gca,'Fontsize',15,'FontWeight','bold')
% % title('Maxiumum Voltage Magnitude Error','Fontsize',15,'FontWeight','bold')
% % xlabel('Substation power |S_{0}|','Fontsize',15,'FontWeight','bold')
% % ylabel('Magnitude error [p.u.]','Fontsize',15,'FontWeight','bold')
% % axis([0 1.5 0 inf])
% % pbaspect([1 0.25 1])
% 
% 
% figure,
% % plot(SS(VV == 0),180/pi*DEM2(VV == 0),'g.',SS(VV == 1),180/pi*DEM2(VV == 1),'r.','LineWidth',2)
% plot(SS,180/pi*DEM2,'r.','LineWidth',2)
% set(gca,'YTick',0:0.1:0.5,'Fontsize',15,'FontWeight','bold')
% title('Maxiumum Voltage Angle Error','Fontsize',15,'FontWeight','bold')
% xlabel('Substation power |S_{0}|','Fontsize',15,'FontWeight','bold')
% ylabel('Angle error [deg]','Fontsize',15,'FontWeight','bold')
% axis([0 1.5 0 0.5])
% pbaspect([1 0.375 1])
% 
% print('-f3','-depsc',['C:\Users\Michael\Desktop\temp\mcangleerror2.eps'])
% print('-f3','-dpng',['C:\Users\Michael\Desktop\temp\mcangleerror2.png'])
% 
% 
% % figure,
% % plot(SS(VV == 0),DEM2P(VV == 0),'g.',SS(VV == 1),DEM2P(VV == 1),'r.','LineWidth',2)
% % set(gca,'Fontsize',15,'FontWeight','bold')
% % title('Maxiumum Voltage Magnitude Error','Fontsize',15,'FontWeight','bold')
% % xlabel('Substation power |S_{0}|','Fontsize',15,'FontWeight','bold')
% % ylabel('Magnitude error [p.u.]','Fontsize',15,'FontWeight','bold')
% % axis([0 1.5 0 inf])
% % pbaspect([1 0.25 1])
% 
% 
% figure,
% % plot(SS,TVEM1,'g.',SS,TVEM2,'r.','LineWidth',2)
% plot(SS,TVEM1,'g.','LineWidth',2)
% set(gca,'Fontsize',15,'FontWeight','bold')
% title('Maxiumum Total Vector Error','Fontsize',15,'FontWeight','bold')
% xlabel('Substation power |S_{0}|','Fontsize',15,'FontWeight','bold')
% ylabel('Vector error [p.u.]','Fontsize',15,'FontWeight','bold')
% axis([0 1.5 0 0.015])
% pbaspect([1 0.375 1])
% 
% print('-f4','-depsc',['C:\Users\Michael\Desktop\temp\mctverror1.eps'])
% print('-f4','-dpng',['C:\Users\Michael\Desktop\temp\mctverror1.png'])
% 
% figure,
% plot(SS,TVEM2,'r.','LineWidth',2)
% set(gca,'Fontsize',15,'FontWeight','bold')
% title('Maxiumum Total Vector Error','Fontsize',15,'FontWeight','bold')
% xlabel('Substation power |S_{0}|','Fontsize',15,'FontWeight','bold')
% ylabel('Vector error [p.u.]','Fontsize',15,'FontWeight','bold')
% axis([0 1.5 0 0.015])
% pbaspect([1 0.375 1])
% 
% print('-f5','-depsc',['C:\Users\Michael\Desktop\temp\mctverror2.eps'])
% print('-f5','-dpng',['C:\Users\Michael\Desktop\temp\mctverror2.png'])
% 
% figure,
% plot(SS,SEM1,'b.','LineWidth',2)
% set(gca,'Fontsize',15,'FontWeight','bold')
% title('Maxiumum Complex Power Error','Fontsize',15,'FontWeight','bold')
% xlabel('Substation power |S_{0}|','Fontsize',15,'FontWeight','bold')
% % ylabel('Power error [p.u.]','Fontsize',15,'FontWeight','bold')
% ylabel([char(949) '_{power} [p.u.]'],'Fontsize',15,'FontWeight','bold')
% axis([0 1.5 0 0.075])
% pbaspect([1 0.375 1])
% 
% print('-f6','-depsc',['C:\Users\Michael\Desktop\temp\mcpowererror.eps'])
% print('-f6','-dpng',['C:\Users\Michael\Desktop\temp\mcpowererror.png'])


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% ylabel([char(949) '_{power} [p.u.]'],'Fontsize',15,'FontWeight','bold')


figure,
plot(SS,EEM1,'r.','LineWidth',2)
set(gca,'XTick',0:0.25:1.5,'YTick',0:0.002:0.02,'Fontsize',12,'FontWeight','bold')
title('Maximum Network Voltage Magnitude Error','Fontsize',12,'FontWeight','bold')
% xlabel('Substation power \Sigma_{\phi \in \{a,b,c\}} |S_{\infty, 650}^{\phi}| [p.u.]','Fontsize',12,'FontWeight','bold')
xlabel('Substation power S_{sub} [p.u.]','Fontsize',12,'FontWeight','bold')
% ylabel('Magnitude error [p.u.]','Fontsize',12,'FontWeight','bold')
ylabel([char(949) '_{mag} [p.u.]'],'Fontsize',12,'FontWeight','bold')
axis([0 1.5 0 0.01])
pbaspect([1 0.25 1])

print('-f1','-depsc',['C:\Users\Michael\Desktop\temp\montecarlomagnitude.eps'])

figure,
plot(SS,180/pi*DEM1,'g.','LineWidth',2)
set(gca,'XTick',0:0.25:1.5,'YTick',0:0.1:0.5,'Fontsize',12,'FontWeight','bold')
title('Maximum Network Voltage Angle Error','Fontsize',12,'FontWeight','bold')
% xlabel('Substation power \Sigma_{\phi \in \{a,b,c\}} |S_{\infty, 650}^{\phi}| [p.u.]','Fontsize',12,'FontWeight','bold')
xlabel('Substation power S_{sub} [p.u.]','Fontsize',12,'FontWeight','bold')
% ylabel('Magnitude error [p.u.]','Fontsize',12,'FontWeight','bold')
ylabel([char(949) '_{angle} [' char(176) ']'],'Fontsize',12,'FontWeight','bold')
axis([0 1.5 0 0.5])
pbaspect([1 0.25 1])

print('-f2','-depsc',['C:\Users\Michael\Desktop\temp\montecarloangle.eps'])

figure,
plot(SS,SEM1,'b.','LineWidth',2)
set(gca,'XTick',0:0.25:1.5,'YTick',0:0.02:0.1,'Fontsize',12,'FontWeight','bold')
title('Maximum Network Apparent Power Error','Fontsize',12,'FontWeight','bold')
% xlabel('Substation power \Sigma_{\phi \in \{a,b,c\}} |S_{\infty, 650}^{\phi}| [p.u.]','Fontsize',12,'FontWeight','bold')
xlabel('Substation power S_{sub} [p.u.]','Fontsize',12,'FontWeight','bold')
% ylabel('Magnitude error [p.u.]','Fontsize',12,'FontWeight','bold')
ylabel([char(949) '_{power} [p.u.]'],'Fontsize',12,'FontWeight','bold')
axis([0 1.5 0 0.1])


pbaspect([1 0.25 1])

print('-f3','-depsc',['C:\Users\Michael\Desktop\temp\montecarlopower.eps'])