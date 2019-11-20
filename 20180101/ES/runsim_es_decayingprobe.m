% Michael Sankur - msankur@lbl.gov
% 2016.08.01

clc, clear all, close all

% path(path,genpath('C:\Users\Michael\Desktop\reconfiguration\Mesh\GEN'));
path(path,genpath('/home/michael/Dropbox/Unbalanced LinDistflow/20180101/'));


% if exist('mxlpsolve','file') == 0
%     path(path,'C:\Users\Michael\Desktop\mxlp');
% end
% 
% if exist('cvx_begin','file') == 0
%     cd C:\Users\Michael\Desktop\cvx-w64\cvx
%     cvx_setup
% end

%% Load feeder

networkname = '05node_singlephase_radial';

fp = 'C:\Users\Michael\Desktop\reconfiguration\Feeders\';
fp = '/home/michael/Dropbox/Unbalanced LinDistflow/20180101/NETWORKS/';
fn = [networkname '.txt'];

[network1] = network_mapper_function_20180101(networkname, fp, fn);

%% Network paramaters

nnode1 = network1.nodes.nnode;
nline1 = network1.lines.nline;

%% Load parameters

network1.loads.aPQ = 0.85*ones(3,nnode1).*network1.nodes.PH;
network1.loads.aI = zeros(3,nnode1);
network1.loads.aZ = 0.15*ones(3,nnode1).*network1.nodes.PH;

% network1.loads.aPQ = 1.00*ones(3,nnode1).*network1.nodes.PH;
% network1.loads.aI = zeros(3,nnode1);
% network1.loads.aZ = 0.00*ones(3,nnode1).*network1.nodes.PH;

network1.loads.spu = 1.5*network1.loads.spu;
network1.loads.ppu = real(network1.loads.spu);
network1.loads.qpu = imag(network1.loads.spu);

%% Capacitor parameters

network1.caps.cappu = 1*network1.caps.cappu;

%% Controller parameters

network1.cons.wmaxpu = 0.5*network1.cons.wmaxpu;

%% Simulation parameters

Vnom = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim.Vnom = Vnom;

sim.VNR = [];
sim.Lmn = [];
sim.Hmn = [];

sim.rho = 0.1;

%%

network1.cons.wpu = zeros(3,nnode1);
[VNR0, INR0, STXNR0, SRXNR0, iNR0, sNR0, iter0] = NR3(network1,[],[],1,Vnom);

%%

dt = 0.01;
time = 0:dt:200;

% vectors
MD = zeros(3,length(time));
ED = zeros(3,length(time));
AD = zeros(3,length(time));
PD = zeros(3,length(time));

P0 = zeros(1,length(time));

J = zeros(1,length(time)); % objective function
Javg = zeros(3,nnode1,length(time));
Jdot = zeros(1,length(time)); % time derivative of objective function
rho = zeros(3,nnode1,2,length(time)); % value after highpass (washout) filter
rhoavg = zeros(3,nnode1,2,length(time));
sigma = zeros(3,nnode1,2,length(time)); % value after demodulation
sigmaavg = zeros(3,nnode1,2,length(time));
xi = zeros(3,nnode1,2,length(time)); % value after lowpass filter
xiavg = zeros(3,nnode1,2,length(time));
uhat = zeros(3,nnode1,length(time)); % value after integration
uhatdot = zeros(3,nnode1,length(time)); % value after integration
u = zeros(3,nnode1,length(time));  % control outputs of ES control (position)
vhat = zeros(3,nnode1,length(time)); % value after integration
vhatdot = zeros(3,nnode1,length(time)); % value after integration
v = zeros(3,nnode1,length(time));
a = zeros(3,nnode1,2,length(time));  % adaptive probing amplitude
sw = zeros(3,nnode1,2,length(time));
% e = zeros(2,length(time),nc);  % LPF of objective


% initial conditions
u(:,:,1) = 0;
v(:,:,1) = 0;
MD(:,1) = abs(VNR0(:,4)) - abs(VNR0(:,6));
ED(:,1) = abs(VNR0(:,4)).^2 - abs(VNR0(:,6)).^2;
AD(:,1) = angle(VNR0(:,4)) - angle(VNR0(:,6));
PD(:,1) = abs(VNR0(:,4)- VNR0(:,6));

% system state vectors
ii = zeros(3,nnode1,length(time));
ss = zeros(3,nnode1,length(time));

%

Ptrack = zeros(1,length(time));
% Ptrack(1,1:60001) = 0.36;
% Ptrack(1,60002:120001) = 0.36 + (0.40 - 0.36)/(time(120001) - time(60002))*(time(60002:120001) - time(60002));
% Ptrack(1,120002:180001) = 0.40;
Ptrack(1,1:10001) = 0.40;
Ptrack(1,10002:20001) = 0.34;

%
J(1) = (real(STXNR0(1,1)) - Ptrack(1))^2;

% network states
VV = VNR0;
II = INR0;
STXSTX = STXNR0;
SRXSRX = SRXNR0;

% setup parameters for the ES controller
a0 = 0.005;
wmaxpu = network1.cons.wmaxpu;
network1.cons.fes(1,:) = [0 0 1 1+(sqrt(2)-1)/10 0 1+(sqrt(3)-1)/10];
fes = network1.cons.fes
hpfes = network1.cons.hpfes
lpfes = network1.cons.lpfes
network1.cons.kintes = 2.5*network1.cons.kintes;
kintes = network1.cons.kintes

betau = zeros(3,nnode1);
betau(network1.cons.fes ~= 0) = 0.1;

betav = zeros(3,nnode1);
betav(network1.cons.fes ~= 0) = 0.1;

%
slacknode = 1;
Vnom = 1*[1*exp(0*1j*pi/180);
    1*exp(-120*1j*pi/180);
    1*exp(120*1j*pi/180)];

%%

tic

for kt = 1:length(time)
    
    % display time
    if mod(kt-1,1000) == 0
        disp(['kt = ' num2str(kt) ' | time = ' num2str(time(kt))])
    end
    
    % start ES control on second timestep
    if kt >= 2
        
        % control from previous timestep
        wk = 1*u(:,:,kt-1) + 1j*v(:,:,kt-1);
        network1.cons.wpu = wk;
        
        % solve power flow
        [Vk, Ik, STXk, SRXk, ik, sk, iterk] = NR3(network1,[],[],slacknode,Vnom);
        
        % network state
        VV(:,:,end+1) = Vk;
        II(:,:,end+1) = Ik;
        STXSTX(:,:,end+1) = STXk;
        SRXSRX(:,:,end+1) = SRXk;
        
        % objective function
        Jk = 1*(real(STXk(1,1)) - Ptrack(kt))^2 + 0*imag(STXk(1,1))^2 + 0*sum(sum(v(:,:,kt-1).^2));
        J(kt) = Jk;
                
        % obtain control
        for ph = 1:3
            for kn = 1:nnode1
                if fes(ph,kn) ~= 0

                    % real power ES loop
                    [uk, rhok, sigmak, xik, uhatk, ak, ~, uhatdotk] = ...
                        esfunction(time(kt), dt, hpfes(ph,kn), lpfes(ph,kn), fes(ph,kn), ...
                        kintes(ph,kn), a(ph,kn,1,kt-1), J(kt), J(kt-1), Jdot(kt-1), ...
                        rho(ph,kn,1,kt-1), sigma(ph,kn,1,kt-1), xi(ph,kn,1,kt-1), uhat(ph,kn,kt-1), a0, 0, sw(ph,kn,1,kt-1), betau(ph,kn));

                    % store internal states
                    rho(ph,kn,1,kt) = rhok;
                    temprho(:,1) = rho(ph,kn,1,:);
                    rhoavg(ph,kn,1,kt) = avgfunction(kt, time(kt), dt, fes(ph,kn), temprho);
                    
                    sigma(ph,kn,1,kt) = sigmak;
                    tempsigma(:,1) = sigma(ph,kn,1,:);
                    sigmaavg(ph,kn,1,kt) = avgfunction(kt, time(kt), dt, fes(ph,kn), tempsigma);
                    
                    xi(ph,kn,1,kt) = xik;
                    tempxi(:,1) = xi(ph,kn,1,:);
                    xiavg(ph,kn,1,kt) = avgfunction(kt, time(kt), dt, fes(ph,kn), tempxi);
                    
                    a(ph,kn,1,kt) = ak;
                    
                    uhatdot(ph,kn,kt) = uhatdotk;

                    % reactive power ES loop
                    [vk, rhok, sigmak, xik, vhatk, ak ~, vhatdotk] = ...
                        esfunction(time(kt), dt, hpfes(ph,kn), lpfes(ph,kn), fes(ph,kn), ...
                        kintes(ph,kn), a(ph,kn,2,kt-1), J(kt), J(kt-1), Jdot(kt-1), ...
                        rho(ph,kn,2,kt-1), sigma(ph,kn,2,kt-1), xi(ph,kn,2,kt-1), vhat(ph,kn,kt-1), a0, -pi/2, sw(ph,kn,2,kt-1),betav(ph,kn));

                    % store internal states
                    rho(ph,kn,2,kt) = rhok;
                    temprho(:,1) = rho(ph,kn,2,:);
                    rhoavg(ph,kn,2,kt) = avgfunction(kt, time(kt), dt, fes(ph,kn), temprho);
                    
                    sigma(ph,kn,2,kt) = sigmak;
                    tempsigma(:,1) = sigma(ph,kn,2,:);
                    sigmaavg(ph,kn,2,kt) = avgfunction(kt, time(kt), dt, fes(ph,kn), tempsigma);
                    
                    xi(ph,kn,2,kt) = xik;
                    tempxi(:,1) = xi(ph,kn,2,:);
                    xiavg(ph,kn,2,kt) = avgfunction(kt, time(kt), dt, fes(ph,kn), tempxi);
                    
                    a(ph,kn,2,kt) = ak;
                    
                    vhatdot(ph,kn,kt) = vhatdotk;

                    % rectify average conrol (whatk) if control (wk) magnitude
                    % is greater than controller apparent power capacity
                    [uhatk,vhatk] = rectify_what_wmax(time(kt),uk,vk,uhatk,vhatk,fes(ph,kn),a0,wmaxpu(ph,kn));

                    % rectify average conrol (whatk) if control (wk) has
                    % positive real power component (consuming real power)
                    [uhatk,vhatk] = rectify_what_generation(time(kt),uk,vk,uhatk,vhatk,fes(ph,kn),a0);
                    
                    uhat(ph,kn,kt) = uhatk;
                    vhat(ph,kn,kt) = vhatk;

                    % recompute current real and reactive power control and
                    % store
                    uk = uhatk + a(ph,kn,1,kt)*cos(2*pi*fes(ph,kn)*time(kt));
                    u(ph,kn,kt) = uk;
                    
                    vk = vhatk + a(ph,kn,2,kt)*sin(2*pi*fes(ph,kn)*time(kt));
                    v(ph,kn,kt) = vk;
                    
                    if kt*dt > 1/min(min(fes(fes ~= 0)))
                        if abs(xiavg(ph,kn,1,kt)/a(ph,kn,1,kt)) <= 0.005 && abs(uhatdotk) <= 5e-5
                            sw1 = 1;
                            sw(ph,kn,1,kt) = 1;
                        elseif Jk >= 3e-4
                            sw(ph,kn,1,kt) = 0;
                        end

                        if abs(xiavg(ph,kn,2,kt)/a(ph,kn,2,kt)) <= 0.005 && abs(vhatdotk) <= 5e-5
                            sw2 = 1;
                            sw(ph,kn,2,kt) = 1;
                        elseif Jk >= 3e-4
                            sw(ph,kn,2,kt) = 0;
                        end
                    end
                
                end                                
            end
        end
    end
end

toc

%%

clear km kn k1 k2 k3 k4 k5 k6 k7 k8 k9
close all

figure, box on, hold on
plot(time,J,'r','LineWidth',1.5)
% Javgtemp(:,1) = Javg(1,3,:);
% plot(time,Javgtemp,'g--')
% Javgtemp(:,1) = Javg(1,4,:);
% plot(time,Javgtemp,'b--')
% Javgtemp(:,1) = Javg(1,6,:);
% plot(time,Javgtemp,'k--')
set(gca,'FontSize',12,'FontWeight','bold','XTick',0:time(end)/4:time(end))
axis([time(1) time(end) 0 inf])

figure, box on, hold on
PTXplot(:,1) = real(STXSTX(1,1,:));
QTXplot(:,1) = imag(STXSTX(1,1,:));
STXplot(:,1) = abs(STXSTX(1,1,:));
plot(time,PTXplot,'r','LineWidth',1.5)
plot(time,Ptrack,'r--','LineWidth',1.5)
plot(time,QTXplot,'g','LineWidth',1.5)
plot(time,STXplot,'b','LineWidth',1.5)
set(gca,'FontSize',12,'FontWeight','bold','XTick',0:time(end)/4:time(end))
axis([time(1) time(end) 0 0.75])

%%

% abcvec = {'a','b','c'};
% for kn = 1:nnode1
%     for ph = 1:3    
%         if network1.cons.wmaxpu(ph,kn) ~= 0
%             figure, box on, hold on
%             
%             subplot(3,1,1), hold on
%             utemp(1,:) = u(ph,kn,:);
%             plot(time,utemp,'r-','LineWidth',1.5)
%             vtemp(1,:) = v(ph,kn,:);
%             plot(time,vtemp,'g-','LineWidth',1.5)
%             plot(time,sqrt(utemp.^2 + vtemp.^2),'b-','LineWidth',1.5)
%             set(gca,'FontSize',12,'FontWeight','bold','XTick',0:time(end)/4:time(end),'YTick',-0.15:0.05:0.15)
%             title(['Control on Phase ' abcvec{ph} ' at Node ' network1.nodes.nodelist{kn}],'FontSize',12,'FontWeight','bold')
% %             xlabel('Time [s]','FontSize',12,'FontWeight','bold')
%             ylabel('[p.u.]','FontSize',12,'FontWeight','bold')
%             legend({'u_{m}(t)','v_{m}(t)','|w_{m}(t)|'},'FontWeight','bold','location','southeast')
%             axis([time(1) time(end) -0.06 0.06])
% %             pbaspect([1 1.25 1])
% %             print(['-f' num2str(5+kn)],'-depsc',['C:\Users\Michael\Desktop\temp\es\' 'w' nodes.nodelist{connode(kn)} '.eps'])
% 
%             subplot(3,1,2), hold on
%             sigmaavgplot1(:,1) = sigmaavg(1,kn,1,:);
%             sigmaavgplot2(:,1) = sigmaavg(1,kn,2,:);            
%             plot(time,sigmaavgplot1,'r','LineWidth',1.5)
%             plot(time,sigmaavgplot2,'g','LineWidth',1.5)
%             set(gca,'FontSize',12,'FontWeight','bold','XTick',0:time(end)/4:time(end))
%             
%             subplot(3,1,3), hold on
%             xiavgplot1(:,1) = xiavg(1,kn,1,:);
%             xiavgplot2(:,1) = xiavg(1,kn,2,:);
%             plot(time,xiavgplot1,'r','LineWidth',1.5)
%             plot(time,xiavgplot2,'g','LineWidth',1.5)
%             set(gca,'FontSize',12,'FontWeight','bold','XTick',0:time(end)/4:time(end))
%            
%         end
%     end
% end

%%

abcvec = {'a','b','c'};
for kn = 1:nnode1
    for ph = 1:3
        if network1.cons.wmaxpu(ph,kn) ~= 0
            figure, box on, hold on
            
            subplot(4,1,1), box on, hold on
            utemp(1,:) = u(ph,kn,:);
            plot(time,utemp,'r-','LineWidth',1.5)
            uhattemp(1,:) = uhat(ph,kn,:);
            plot(time,uhattemp,'r--','LineWidth',1.5)
            vtemp(1,:) = v(ph,kn,:);
            plot(time,vtemp,'g-','LineWidth',1.5)
            vhattemp(1,:) = vhat(ph,kn,:);
            plot(time,vhattemp,'g--','LineWidth',1.5)
            plot(time,sqrt(utemp.^2 + vtemp.^2),'b-','LineWidth',1.5)
            plot(time,sqrt(uhattemp.^2 + vhattemp.^2),'b-','LineWidth',1.5)
            set(gca,'FontSize',12,'FontWeight','bold','XTick',0:time(end)/4:time(end),'YTick',-0.15:0.05:0.15)
            title(['Control on Phase ' abcvec{ph} ' at Node ' network1.nodes.nodelist{kn}],'FontSize',12,'FontWeight','bold')
%             xlabel('Time [s]','FontSize',12,'FontWeight','bold')
            ylabel('[p.u.]','FontSize',12,'FontWeight','bold')
%             legend({'u_{m}(t)','v_{m}(t)','|w_{m}(t)|'},'FontWeight','bold','location','southeast')
            axis([time(1) time(end) -0.06 0.06])
%             pbaspect([1 1.25 1])
%             print(['-f' num2str(5+kn)],'-depsc',['C:\Users\Michael\Desktop\temp\es\' 'w' nodes.nodelist{connode(kn)} '.eps'])

%             subplot(4,1,2), hold on
%             sigmaavgplot1(:,1) = sigmaavg(1,kn,1,:)./a(1,kn,1,:);
%             sigmaavgplot2(:,1) = sigmaavg(1,kn,2,:)./a(1,kn,2,:);            
%             plot(time,sigmaavgplot1,'r','LineWidth',1.5)
%             plot(time,sigmaavgplot2,'g','LineWidth',1.5)
%             set(gca,'FontSize',12,'FontWeight','bold','XTick',0:time(end)/4:time(end))

            subplot(4,1,2), box on, hold on
            xiavgplot1(:,1) = xiavg(1,kn,1,:)./a(1,kn,1,:);
            xiavgplot2(:,1) = xiavg(1,kn,2,:)./a(1,kn,2,:);
            plot(time,xiavgplot1,'r','LineWidth',1.5)
            plot(time,xiavgplot2,'g','LineWidth',1.5)
            plot(time,zeros(size(time)),'k--','LineWidth',1.5)
            set(gca,'FontSize',12,'FontWeight','bold','XTick',0:time(end)/4:time(end))
            title('\xi_{avg}','FontSize',12,'FontWeight','bold','Interpreter','latex')

            subplot(4,1,3), box on, hold on
            uhatdotplot(:,1) = uhatdot(ph,kn,:);
            vhatdotplot(:,1) = vhatdot(ph,kn,:);
            plot(time,uhatdotplot,'r','LineWidth',1.5)
            plot(time,vhatdotplot,'g','LineWidth',1.5)
            plot(time,zeros(size(time)),'k--','LineWidth',1.5)
            set(gca,'FontSize',12,'FontWeight','bold','XTick',0:time(end)/4:time(end))    
            title('uhatdot or vhatdot','FontSize',12,'FontWeight','bold')%,'Interpreter','latex')
                        
            subplot(4,1,4), box on, hold on
            swplot1(:,1) = sw(1,kn,1,:);
            swplot2(:,1) = sw(1,kn,2,:);
            plot(time,swplot1,'r-','LineWidth',1.5)
            plot(time,swplot2,'g--','LineWidth',1.5)
            set(gca,'FontSize',12,'FontWeight','bold','XTick',0:time(end)/4:time(end))
            title('Switch','FontSize',12,'FontWeight','bold')
           
        end
    end
end

%%

% figure, box on, hold on
% for k1 = 1:nnode1
%    sigmaavgplot1(:,1) = sigmaavg(1,k1,1,:);
%    sigmaavgplot2(:,1) = sigmaavg(1,k1,2,:);
%    plot(time,sigmaavgplot1,'r')
%    plot(time,sigmaavgplot2,'g')
% end
% 
% figure, box on, hold on
% for k1 = 1:nnode1
%    xiavgplot1(:,1) = xiavg(1,k1,1,:);
%    xiavgplot2(:,1) = xiavg(1,k1,2,:);
%    plot(time,xiavgplot1,'r')
%    plot(time,xiavgplot2,'g')
% end
% 
% %%
% 
% abcvec = {'a','b','c'};
% for kn = 1:nnode1
%     for ph = 1:3    
%         if network1.cons.wmaxpu(ph,kn) ~= 0
%             figure, box on, hold on
%             utemp(1,:) = u(ph,kn,:);
%             plot(time,utemp,'r-','LineWidth',1.5)
%             vtemp(1,:) = v(ph,kn,:);
%             plot(time,vtemp,'g-','LineWidth',1.5)
%             plot(time,sqrt(utemp.^2 + vtemp.^2),'b-','LineWidth',1.5)
%             set(gca,'FontSize',12,'FontWeight','bold','XTick',0:time(end)/4:time(end),'YTick',-0.15:0.05:0.15)
%             title(['Control on Phase ' abcvec{ph} ' at Node ' network1.nodes.nodelist{kn}],'FontSize',12,'FontWeight','bold')
%             xlabel('Time [s]','FontSize',12,'FontWeight','bold')
%             ylabel('[p.u.]','FontSize',12,'FontWeight','bold')
%             legend({'u_{m}(t)','v_{m}(t)','|w_{m}(t)|'},'FontWeight','bold','location','southeast')
%             axis([time(1) time(end) -0.15 0.15])
%             pbaspect([1 1.25 1])
% %             print(['-f' num2str(5+kn)],'-depsc',['C:\Users\Michael\Desktop\temp\es\' 'w' nodes.nodelist{connode(kn)} '.eps'])
%         end
%     end
% end
% 
%%

x = 0:0.001:50;
figure, box on, hold on,
y1 = 0;
y2 = 0;
for k1 = 1:nnode1
    if fes(1,k1) ~= 0
        y1 = y1 + cos(2*pi*fes(1,k1)*x);
        y2 = y2 + sin(2*pi*fes(1,k1)*x);
    end
end
plot(x,y1,'r');
plot(x,y2,'g');
plot(x,y1+y2,'b');

