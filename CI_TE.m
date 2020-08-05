clc, clear, close all

%%

A0 = 1;
B0 = 0;

A0 = -1/2
B0 = -sqrt(3)/2

dA = linspace(-0.2,0.2,101);
dB = linspace(-0.2,0.2,101);

f0 = zeros(length(dA),length(dB));
f1 = zeros(length(dA),length(dB));

gf = [(A0/(A0^2 + B0^2)^(1/2)), (B0/(A0^2 + B0^2)^(1/2))];

Hf(1,1) = ((A0^2 + B0^2)^(1/2) - A0^2*(A0^2 + B0^2)^(-1/2))/(A0^2 + B0^2);
Hf(1,2) = (-A0*B0*(A0^2 + B0^2)^(-1/2))/(A0^2 + B0^2);
Hf(2,1) = (-A0*B0*(A0^2 + B0^2)^(-1/2))/(A0^2 + B0^2);
Hf(2,2) = ((A0^2 + B0^2)^(1/2) - B0^2*(A0^2 + B0^2)^(-1/2))/(A0^2 + B0^2);

for k1 = 1:length(dA)
    for k2 = 1:length(dB)
        
        f0(k1,k2) = ((A0+dA(k1))^2 + (B0+dB(k2))^2)^(1/2);
                
        f1(k1,k2) = (A0^2 + B0^2)^(1/2) + gf*[dA(k1); dB(k2)];
        
        f2(k1,k2) = (A0^2 + B0^2)^(1/2) + gf*[dA(k1); dB(k2)] + (1/2)*[dA(k1) dB(k2)]*Hf*[dA(k1); dB(k2)];
        
    end
end
        

figure, box on, hold on
surface(A0+dA,B0+dB,f0,'linestyle','none')
surface(A0+dA,B0+dB,f1,'linestyle','none')
% surface(A0+dA,B0+dB,f2,'linestyle','none')
xlabel('A')
ylabel('B')

max(max(abs(f1-f0)))
max(max(abs(f2-f0)))

%% Taylor Expansion for Voltage Magnitude Term

Vmagnom = 1;
Vangnom = 120;

dVmag = linspace(-0.05,+0.05,101);
dVang = linspace(-1,1,101);

Vmag = Vmagnom + dVmag;
Vang = Vangnom + dVang;

f0 = zeros(length(Vmag),length(Vang));
f1 = zeros(length(Vmag),length(Vang));
f2 = zeros(length(Vmag),length(Vang));

A0 = Vmagnom*cos(pi/180*Vangnom);
B0 = Vmagnom*sin(pi/180*Vangnom);

gf = [(A0/(A0^2 + B0^2)^(1/2)), (B0/(A0^2 + B0^2)^(1/2))];

Hf = zeros(2,2);
Hf(1,1) = ((A0^2 + B0^2)^(1/2) - A0^2*(A0^2 + B0^2)^(-1/2))/(A0^2 + B0^2);
Hf(1,2) = (-A0*B0*(A0^2 + B0^2)^(-1/2))/(A0^2 + B0^2);
Hf(2,1) = (-A0*B0*(A0^2 + B0^2)^(-1/2))/(A0^2 + B0^2);
Hf(2,2) = ((A0^2 + B0^2)^(1/2) - B0^2*(A0^2 + B0^2)^(-1/2))/(A0^2 + B0^2);

for k1 = 1:length(Vmag)
    for k2 = 1:length(Vang)
        
        A = Vmag(k1)*cos(pi/180*Vang(k2));
        B = Vmag(k1)*sin(pi/180*Vang(k2));
        
        dA = A - A0;
        dB = B - B0;
        
        f0(k1,k2) = (A^2 + B^2)^(1/2);
        
        f1(k1,k2) = (A0^2 + B0^2)^(1/2) + gf*[dA; dB];
        
        f2(k1,k2) = (A0^2 + B0^2)^(1/2) + gf*[dA; dB] + (1/2)*[dA, dB]*Hf*[dA; dB];
        
        
    end
end

figure, box on, hold on
surface(Vmag,Vang,f0,'linestyle','none')
surface(Vmag,Vang,f1,'linestyle','none')
surface(A0+dA,B0+dB,f2,'linestyle','none')
xlabel('A')
ylabel('B')

disp('Max Error First Order Magnitude')
max(max(abs(f1-f0)))
disp('Max Error Second Order Magnitude')
max(max(abs(f2-f0)))

%% Taylor Expansion for Squared Voltage Magnitude Term

Vmagnom = 1;
Vangnom = 120;

dVmag = linspace(-0.05,+0.05,101);
dVang = linspace(-1,1,101);

Vmag = Vmagnom + dVmag;
Vang = Vangnom + dVang;

f0 = zeros(length(Vmag),length(Vang));
f1 = zeros(length(Vmag),length(Vang));
f2 = zeros(length(Vmag),length(Vang));

A0 = Vmagnom*cos(pi/180*Vangnom);
B0 = Vmagnom*sin(pi/180*Vangnom);

gf = [2*A0, 2*B0];

Hf = [2, 0;
    0, 2];

for k1 = 1:length(Vmag)
    for k2 = 1:length(Vang)
        
        A = Vmag(k1)*cos(pi/180*Vang(k2));
        B = Vmag(k1)*sin(pi/180*Vang(k2));
        
        dA = A - A0;
        dB = B - B0;
        
        f0(k1,k2) = (A^2 + B^2);
        
        f1(k1,k2) = (A0^2 + B0^2) + gf*[dA; dB];
        
        f2(k1,k2) = (A0^2 + B0^2)^(1/2) + gf*[dA; dB] + (1/2)*[dA, dB]*Hf*[dA; dB];
        
        
    end
end

figure, box on, hold on
surface(Vmag,Vang,f0,'linestyle','none')
surface(Vmag,Vang,f1,'linestyle','none')
surface(A0+dA,B0+dB,f2,'linestyle','none')
xlabel('A')
ylabel('B')

disp('Max Error First Order Squared Magnitude')
max(max(abs(f1-f0)))
disp('Max Error Second Order Squared Magnitude')
max(max(abs(f2-f0)))
