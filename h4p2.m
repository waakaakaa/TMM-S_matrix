clear
clc

degrees = pi/180;

nPoints = 1000;
lam0 = linspace(1e-6,2e-6,nPoints); %free space wavelength
Rvector = zeros(1,nPoints);
Tvector = zeros(1,nPoints);

theta = 0 * degrees; %elevation angle
phi = 0 * degrees; %azimuthal angle
pte = 1; %amplitude of TE polarization
ptm = 0; %amplitude of TM polarization

ur1 = 1; %permeability in the reflection region
er1 = 1; %permittivity in the reflection region
ur2 = 1; %permeability in the transmission region
er2 = 1; %permittivity in the transmission region
UR = [ 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00]; %array of permeabilities in each layer
ER = [ 2.25 4.41 2.25 4.41 2.25 4.41 2.25 4.41 2.25 4.41 2.25 4.41 2.25 4.41 2.25 4.41 2.25 4.41 2.25 4.41]; %array of permittivities in each layer
L =  [ 250  180  250  180  250  180  250  180  250  180  250  180  250  180  250  180  250  180  250  180]*1e-9; %array of the thickness of each layer

DEV = {er1,ur1,er2,ur2,ER,UR,L};

count = 1;
for i = lam0
    SRC = {i,theta,phi,pte,ptm};
    DAT = tmm1d(DEV,SRC);
    Rvector(1,count) = DAT{1};
    Tvector(1,count) = DAT{2};
    count = count + 1;
end
subplot(2,1,1);hold on;box on;
plot(lam0*1e6,Rvector,'r');
plot(lam0*1e6,Tvector,'b--');
subplot(2,1,2);hold on;box on;
plot(lam0*1e6,10*log10(Rvector),'r');
plot(lam0*1e6,10*log10(Tvector),'b--');
ylim([-30 0])