% Michael Sankur - msankur@berkeley.edu
% 2016.06.03

% This is a one time step simulation in which DERs are dispatched to
% balance the voltage across the existing phases on each node of the
% feeder.

clc, clear all, close all

if exist('FBS_3phase_fun_20160603','file') == 0 || exist('feeder_mapper_ieee_function_20160603','file') == 0
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

feedername = 'ieee_13node_balance';

fp = 'C:\Users\Michael\Desktop\reconfiguration\Feeders\';
fn = [feedername '.txt'];

[feeder, nodes, lines, configs, loads, caps, controllers] = feeder_mapper_TX_function_20170215(feedername, fp, fn);

%% Feeder paramaters

radflag1 = 1;
for k1 = 1:size(nodes.FM,1)
    if sum(nodes.FM(k1,:) == -1) > 1
        radflag1 = 0;
    end    
end
radflag1

nnode = nodes.nnode;
nline = lines.nline;

%% Load parameters

loads.aPQ = 0.85*ones(3,nnode).*nodes.PH;
loads.aI = zeros(3,nnode);
loads.aZ = 0.15*ones(3,nnode).*nodes.PH;

% loads.spu = 1.125*loads.spu;

%% Capacitor parameters

caps.cappu = 0*caps.cappu;

%% Controller parameters

controllers.wmaxpu = 0.25*controllers.wmaxpu;

%% Simulation parameters

Vnom = 1*[1;
    1*exp(j*-120*pi/180);
    1*exp(j*120*pi/180)];
sim.Vnom = Vnom;

sim.Vfbs = [];
sim.Lfbs = [];
sim.Hfbs = [];

slacknode = 1;
sim.slacknode = 1;

sim.rho = 0.1;

%%

ngrid = 16;

Vmagmin = zeros(3,nnode); Vmagmin(nodes.PH == 0 ) = NaN;
Vmagmax = zeros(3,nnode); Vmagmax(nodes.PH == 0 ) = NaN;
Vmaggrid = zeros(3,ngrid,nnode);
thetamin = zeros(3,ngrid,nnode);
thetamax = zeros(3,ngrid,nnode);

%%

[nvar, Aineq, bineq, Aeq, beq, Aw, bw] = createCVXMatrixFeasible(feeder, nodes, lines, configs, loads, caps, controllers, sim);

%% Find max and min y values

disp('COMPUTE MIN AND MAX MAGNITUDE')

ngrid = 11;

sim.tph = 1;
sim.tnode = 11;

sim.gp = [-1 -1 -1];

[optsolmin, optsolmax] = Solver_LinDist3Flow_TX_Feasible_Eminmax(feeder, nodes, lines, configs, loads, caps, controllers, sim);

agrid = linspace(optsolmin.Vmin,optsolmax.Vmax,ngrid+1);

bmin = [];
bmax = [];
cmin = [];
cmax = [];

for ka = 1:length(agrid)
    
    sim.tph = 2;
    sim.gp = [agrid(ka) -1 -1];

    [optsolmin, optsolmax] = Solver_LinDist3Flow_TX_Feasible_Eminmax(feeder, nodes, lines, configs, loads, caps, controllers, sim);
    
    [optsolmin.Vmin optsolmax.Vmax];
    
    bmin(ka) = optsolmin.Vmin;
    bmax(ka) = optsolmax.Vmax;
    
    bgrid = linspace(optsolmin.Vmin,optsolmax.Vmax,ngrid+1);
    
    for kb = 1:length(bgrid)
        
        sim.tph = 3;
        sim.gp = [agrid(ka) bgrid(kb) -1];
        
        [optsolmin, optsolmax] = Solver_LinDist3Flow_TX_Feasible_Eminmax(feeder, nodes, lines, configs, loads, caps, controllers, sim);
        
        cmin(ka,kb) = optsolmin.Vmin;
        
        cmax(ka,kb) = optsolmax.Vmax;
        
    end
    
end

%%

close all

figure, box on, hold on
plot(agrid,bmin,'b.',agrid,bmax,'g.','Linewidth',2)


% %% Find max and min theta values
% 
% disp('COMPUTE MIN AND MAX THETA')
% 
% sim.Vmaggrid = Vmaggrid;
% 
% % for knode = 1:n
%     for ph = 1:3
%         tic
%         if feeder.PH(ph,knode) == 1
%             
%             disp(['Node: ' num2str(knode) ' - ' char(feeder.nodelist(knode))])
%             disp(['Phase: ' num2str(ph)])
%             
%             sim.tn = knode;
%             sim.tp = ph; 
%             
%             for kgrid = 1:ngrid
%                 sim.iter = 0;
% 
%                 sim.Vfbs = [];
%                 sim.Lfbs = [];
%                 sim.Hfbs = [];
%                 
%                 sim.kgrid = kgrid;
% 
%                 [OPTSOLTHETAMIN, OPTSOLTHETAMAX] = Solver_LinDist3Flow_thetaminmax_20160603(feeder,loads,controllers,sim);
%                 
%                 thetamin(ph,kgrid,knode) = OPTSOLTHETAMIN.Ztheta;
%                 
%                 thetamax(ph,kgrid,knode) = OPTSOLTHETAMAX.Ztheta;
% 
%             end
%             disp(' ')
%         end
%         toc
%     end
% % end
% 
% %% Plot convex hull of points
% 
% clear X Xa Xb Xc Y Ya Yb Yc
% close all
% 
% abc = {'a','b','c'};
% 
% for ph = 1:3
%     if feeder.PH(ph,knode) == 1
%         Xnode = [Vmaggrid(ph,:,knode) Vmaggrid(ph,:,knode)]';
%         Ynode = [thetamin(ph,:,knode) thetamax(ph,:,knode)]';
%         Knode = convhull(Xnode,Ynode);
%         
%         figure,
%         plot(abs(FBSbase.V(ph,knode)),180/pi*angle(FBSbase.V(ph,knode)),'k*', ...
%             Xnode,Ynode,'b.',Xnode(Knode),Ynode(Knode),'r-','MarkerSize',15,'LineWidth',1.5)
%         set(gca,'XTick',0.95:0.01:1.05,'FontSize',15,'FontWeight','bold')
%         title(['Node ' char(feeder.nodelist(knode)) ' phase ' abc{ph}],'FontSize',15,'FontWeight','bold')
%         xlabel(['Magnitude |V_{' char(feeder.nodelist(knode)) '}^{' abc{ph} '}|'],'FontSize',15,'FontWeight','bold')
%         ylabel(['Angle \angle V_{' char(feeder.nodelist(knode)) '}^{' abc{ph} '}'],'FontSize',15,'FontWeight','bold')
% %         axis([0.95 1 -2.2 -0.2])
% 
%         print(['-f' num2str(ph)],'-dpng',['C:\Users\Michael\Desktop\temp\node' char(feeder.nodelist(knode)) abc{ph} '.png'])
%     end
% end
%    
% Xa = [Vmaggrid(1,:,knode) Vmaggrid(1,:,knode)]';
% Ya = [thetamin(1,:,knode) thetamax(1,:,knode)]';
% Ka = convhull(Xa,Ya);
% 
% Xb = [Vmaggrid(2,:,knode) Vmaggrid(2,:,knode)]';
% Yb = [thetamin(2,:,knode) thetamax(2,:,knode)]';
% Kb = convhull(Xb,Yb);
% 
% Xc = [Vmaggrid(3,:,knode) Vmaggrid(3,:,knode)]';
% Yc = [thetamin(3,:,knode) thetamax(3,:,knode)]';
% Kc = convhull(Xc,Yc);
% 
% % figure,
% % plot(abs(FBSbase.V(1,knode)),180/pi*angle(FBSbase.V(1,knode)),'k*', ...
% %     Xa,Ya,'b.',Xa(Ka),Ya(Ka),'r-','MarkerSize',15,'LineWidth',1.5)
% % set(gca,'XTick',0.95:0.01:1.05,'FontSize',15,'FontWeight','bold')
% % title(['Node ' char(feeder.nodelist(knode)) ' phase a'],'FontSize',15,'FontWeight','bold')
% % xlabel(['Magnitude |V_{' char(feeder.nodelist(knode)) '}^{a}|'],'FontSize',15,'FontWeight','bold')
% % ylabel(['Angle \angle V_{' char(feeder.nodelist(knode)) '}^{a}'],'FontSize',15,'FontWeight','bold')
% % 
% % print('-f3','-dpng',['C:\Users\Michael\Desktop\temp\node' char(feeder.nodelist(knode)) 'a.png'])
% % 
% % figure,
% % plot(abs(FBSbase.V(2,knode)),180/pi*angle(FBSbase.V(2,knode)),'k*', ...
% %     Xb,Yb,'b.',Xb(Kb),Yb(Kb),'r-','MarkerSize',15,'LineWidth',1.5)
% % set(gca,'XTick',0.95:0.01:1.05,'FontSize',15,'FontWeight','bold')
% % title(['Node ' char(feeder.nodelist(knode)) ' phase b'],'FontSize',15,'FontWeight','bold')
% % xlabel(['Magnitude |V_{' char(feeder.nodelist(knode)) '}^{b}|'],'FontSize',15,'FontWeight','bold')
% % ylabel(['Angle \angle V_{' char(feeder.nodelist(knode)) '}^{b}'],'FontSize',15,'FontWeight','bold')
% % 
% % print('-f4','-dpng',['C:\Users\Michael\Desktop\temp\node' char(feeder.nodelist(knode)) 'b.png'])
% % 
% % figure,
% % plot(abs(FBSbase.V(3,knode)),180/pi*angle(FBSbase.V(3,knode)),'k*', ...
% %     Xc,Yc,'b.',Xc(Kc),Yc(Kc),'r-','MarkerSize',15,'LineWidth',1.5)
% % set(gca,'XTick',0.95:0.01:1.05,'FontSize',15,'FontWeight','bold')
% % title(['Node ' char(feeder.nodelist(knode)) ' phase c'],'FontSize',15,'FontWeight','bold')
% % xlabel(['Magnitude |V_{' char(feeder.nodelist(knode)) '}^{c}|'],'FontSize',15,'FontWeight','bold')
% % ylabel(['Angle \angle V_{' char(feeder.nodelist(knode)) '}^{c}'],'FontSize',15,'FontWeight','bold')
% % 
% % print('-f5','-dpng',['C:\Users\Michael\Desktop\temp\node'
% % char(feeder.nodelist(knode)) 'c.png'])
% 
% figure, box on, hold on
% plot(Xa(Ka),Ya(Ka),'r.-',Xb(Kb),Yb(Kb)+120,'g.-',Xc(Kc),Yc(Kc)-120,'b.-',...
%     'MarkerSize',15,'LineWidth',1.5)
% plot(abs(FBSbase.V(1,knode)),180/pi*angle(FBSbase.V(1,knode)),'r*',...
%     abs(FBSbase.V(2,knode)),180/pi*angle(FBSbase.V(2,knode))+120,'g*',...
%     abs(FBSbase.V(3,knode)),180/pi*angle(FBSbase.V(3,knode))-120,'b*',...
%     'MarkerSize',15,'LineWidth',1.5)
% set(gca,'XTick',0.90:0.01:1.1,'YTick',-360:0.5:360,'FontSize',15,'FontWeight','bold')
% legend({'a','b','c'},'FontSize',15,'FontWeight','bold')
% 
% print('-f4','-dpng',['C:\Users\Michael\Desktop\temp\' char(feeder.nodelist(knode)) 'feas.png'])
