% This MATLAB program implements the transfer matrix method.
% INITIALIZE MATLAB
close all;
clc;
clear all;
% UNITS
degrees = pi/180;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOURCE PARAMETERS
lam0 = 2.7; %free space wavelength
theta = 57 * degrees; %elevation angle
phi = 23 * degrees; %azimuthal angle
pte = 1/sqrt(2); %amplitude of TE polarization
ptm = 1i/sqrt(2); %amplitude of TM polarization
% EXTERNAL MATERIALS
ur1 = 1.2; %permeability in the reflection region
er1 = 1.4; %permittivity in the reflection region
ur2 = 1.6; %permeability in the transmission region
er2 = 1.8; %permittivity in the transmission region
% DEFINE LAYERS
UR = [ 1 3 ]; %array of permeabilities in each layer
ER = [ 2 1 ]; %array of permittivities in each layer
L = [0.25 0.5]*lam0; %array of the thickness of each layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMPLEMENT TRANSFER MATRIX METHOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('***************************************************************\n');
fprintf('*                             SETUP                           *\n');
fprintf('***************************************************************\n');
k0 = 2*pi/lam0

kx = sqrt(ur1*er1)*sin(theta)*cos(phi)
ky = sqrt(ur1*er1)*sin(theta)*sin(phi)

Qh = [kx*ky 1+ky^2; -(1+kx^2) -kx*ky]
Vh = -1i*Qh

fprintf('***************************************************************\n');
fprintf('*                    ITERATE THROUGH LAYERS                   *\n');
fprintf('***************************************************************\n');

SG11 = zeros(2,2);
SG12 = eye(2,2);
SG21 = eye(2,2);
SG22 = zeros(2,2);
for i = 1:length(L)
    fprintf('-----------------------------------------------> Layer %d of %d\n',i,length(L));
    Q = 1/UR(i) * [kx*ky  UR(i)*ER(i)-kx^2; ky^2-UR(i)*ER(i)  -kx*ky]
    kz = sqrt(UR(i)*ER(i) - kx^2 - ky^2)
    OMEGA = 1i*kz*eye(2,2)
    V = Q * OMEGA^(-1)
    X = expm(OMEGA*k0*L(i))
    A = eye(2,2) + V^(-1)*Vh
    B = eye(2,2) - V^(-1)*Vh
    D = A - X*B*A^(-1)*X*B
    S11 = D^(-1)*(X*B*A^(-1)*X*A - B)
    S12 = D^(-1)*X*(A - B*A^(-1)*B)
    S21 = S12
    S22 = S11
    S = star({SG11 SG12 SG21 SG22},{S11 S12 S21 S22});
    SG11 = S{1}
    SG12 = S{2}
    SG21 = S{3}
    SG22 = S{4}
end
fprintf('***************************************************************\n');
fprintf('*                       EXTERNAL REGIONS                      *\n');
fprintf('***************************************************************\n');
Q = 1/ur1 * [kx*ky  ur1*er1-kx^2; ky^2-ur1*er1  -kx*ky];
kz = sqrt(ur1*er1 - kx^2 - ky^2);
OMEGA = 1i*kz*eye(2,2);
Vref = Q * OMEGA^(-1)
Aref = eye(2,2) + Vh^(-1)*Vref;
Bref = eye(2,2) - Vh^(-1)*Vref;
SR11 = -Aref^(-1)*Bref
SR12 = 2*Aref^(-1)
SR21 = 0.5*(Aref - Bref*Aref^(-1)*Bref)
SR22 = Bref*Aref^(-1)

Q = 1/ur2 * [kx*ky  ur2*er2-kx^2; ky^2-ur2*er2  -kx*ky];
kz = sqrt(ur2*er2 - kx^2 - ky^2);
OMEGA = 1i*kz*eye(2,2);
Vtrn = Q * OMEGA^(-1);
Atrn = eye(2,2) + Vh^(-1)*Vtrn;
Btrn = eye(2,2) - Vh^(-1)*Vtrn;
ST11 = Btrn*Atrn^(-1)
ST12 = 0.5*(Atrn - Btrn*Atrn^(-1)*Btrn)
ST21 = 2*Atrn^(-1)
ST22 = -Atrn^(-1)*Btrn

S = star({SR11 SR12 SR21 SR22},{SG11 SG12 SG21 SG22});
SG11 = S{1};
SG12 = S{2};
SG21 = S{3};
SG22 = S{4};
S = star({SG11 SG12 SG21 SG22},{ST11 ST12 ST21 ST22});
SG11 = S{1}
SG12 = S{2}
SG21 = S{3}
SG22 = S{4}

fprintf('***************************************************************\n');
fprintf('*                   SOLVE SCATTERING PROBLEM                  *\n');
fprintf('***************************************************************\n');
n_inc = sqrt(ur1*er1);
k_inc = k0*n_inc*[sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)]
n_surface_normal = [0;0;1];
if theta == 0
    a_te = [0;1;0]
else
    a_te = cross(n_surface_normal,k_inc) / normest(cross(n_surface_normal,k_inc))
end
a_tm = cross(a_te,k_inc) / normest(cross(a_te,k_inc))
p = pte*a_te + ptm*a_tm

Esrc = [p(1,1);p(2,1)]
Eref = SG11*Esrc
Etrn = SG21*Esrc

Eztrn = - ( kx*Etrn(1,1) + ky*Etrn(2,1) ) /kz
kzref = sqrt(ur1*er1)*cos(theta);
Ezref = - ( kx*Eref(1,1) + ky*Eref(2,1) ) / kzref

R = ( normest([Eref;Ezref]) / normest(p) )^2
T = ( normest([Etrn;Eztrn]) / normest(p) )^2 * real((ur1/ur2)*(kz/kzref))