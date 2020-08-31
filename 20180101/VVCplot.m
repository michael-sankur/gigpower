clc, clear, close all

%%

Vmag = 0.9:0.001:1.1;
E = Vmag.^2;

V1 = 0.95;
V2 = 1.05;

qmin = -0.05;
qmax = 0.05;

for k1 = 1:length(Vmag)
   if Vmag(k1) <= V1
       VVC(k1) = qmin;
   elseif Vmag(k1) >= V1 && Vmag(k1) <= V2
       VVC(k1) = (qmax - qmin)/(V2 - V1)*(Vmag(k1) - V1) + qmin;
   elseif Vmag(k1) >= V2
       VVC(k1) = qmax;
   end
   if E(k1) <= V1^2
       VVC2(k1) = qmin
   elseif E(k1) >= V1^2 && E(k1) <= V2^2
       VVC2(k1) = (qmax - qmin)/(V2 - V1)*((1 + E(k1))/2 - V1) + qmin;
   elseif E(k1) >= V2^2
       VVC2(k1) = qmax;
   end
end

figure, box on, hold on
plot(Vmag,VVC,'b-','LineWidth',2)
% plot(Vmag,VVC2,'g-','LineWidth',2)
set(gca,'XTick',0.9:0.025:1.1,'YTick',qmin:(qmax-qmin)/5:qmax,'Fontsize',12,'FontWeight','bold')
title('Volt VAr Control','Fontsize',12,'FontWeight','bold')
% xlabel('Substation power \Sigma_{\phi \in \{a,b,c\}} |S_{\infty, 650}^{\phi}| [p.u.]','Fontsize',12,'FontWeight','bold')
xlabel('Voltage Magnitude [p.u.]','Fontsize',12,'FontWeight','bold')
ylabel('VAr dispatch [p.u.]','Fontsize',12,'FontWeight','bold')
axis([0.9 1.1 -0.06 0.06])
pbaspect([1 0.25 1])
print('-f1','-depsc',['~/Desktop/temp/VVCplot.eps'])