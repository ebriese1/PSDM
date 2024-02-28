%% Read in data 
clear;clc;
% Read data from excel spreadsheet "PSDM.xlsx"
input = xlsread('PSDM.xlsx');
% Read index No. for multiruns
TestID = input(4,3);

% Read data for model input
NCOMP = input(1,1); % NCOMP: number of components
    % XK = input(6,1:NCOMP); %Freundlich isotherm capacity K,(mg/g)(L/mg)^(1/n)
XK = input(6,1:NCOMP); %Freundlich isotherm capacity K,(mg/g)(L/mg)^(1/n)
XN = input(7,1:NCOMP); % Freundlich isotherm exponents 1/n, dimensionless
Xn = 1./XN; %Freundlich isotherm exponent reciprocal n, dim.
CB0 = input(8,1:NCOMP);% Initial concentrations,ug/L
MW = input(9,1:NCOMP);% Molecular weight, g/mol
VB = input(10,1:NCOMP);% molar volume at the boiling point

L = input(13,1);% L: length of column,cm
DIA = input(14,1); % diameter of column, cm
WT = input(15,1); % WT: mass of adsorbent in column, grams
FLRT = input(16,1); % Flow rate, ml/min
TEMP = input(18,1); % water temperature, celsius degree; 
TEMP = TEMP + 273.15; % water temperature, Kelvin;
TEMPA = TEMP/324.65; % Absolute water temperature, dim.;

RHOP = input(21,1);% apparent adsorbent density, g/cm^3
EPOR = input(22,1);% Particle porosity: fraction of volumetric space in adsorbent phase unoccupied by adsorbent (pore void fraction)
RAD = input(23,1); % Radius of adsorbent, cm
TOR = input(24,1); % tortuosity of particle
% Density of water at a certain temperature, g/cm^3 = 1000 kg/m^3
% Calculated using a correlation developed at Michigan Technological University
RHOL = (-1.41768+8.976651*TEMPA-12.2755*TEMPA^2+7.45844*TEMPA^3-1.738492*TEMPA^4)*0.98396;
% Viscosity of water at a certain temperature,cP (=0.001 kg/m/s)
% Calculated using a correlation shown by Reid (1987) and Yaws (1976)
VISC = exp(-24.71+4209/TEMP+0.04527*TEMP-3.376*10^(-5)*TEMP^2);  

KF = input(26,1:NCOMP); % Film transfer coefficient, cm/s
DS = input(27,1:NCOMP); % Surface diffusivity,cm^2/s
DP = input(28,1:NCOMP);% Pore diffusivity based on pore void fraction,cm^2/s
DL = input(34,1:NCOMP); % Liquid diffusivity, cm^2/s
SC = input(35,1:NCOMP); % Schimit number, dim.
DTOL = input(31,1); % total run time, min.
DSTEP = 1; % Time step added in integration, min.
DT0 = 0.1; % Initial time value in integration, min.

%% Calculate fixed bed parameters
AREA = pi*DIA*DIA/4; % Cross sectional area of differential slice in bed, cm^2
BEDVOL = L*AREA; % Bed volume, ml
EBED = 1.0 - WT/(BEDVOL*RHOP); % Bed porosity: void fraction (Fraction of volumetric space in reactor unoccupied by adsorbent)
EBED = 0.45; 
EBCT = BEDVOL/FLRT;  % Empty bed contact time, min
TAU = EBCT*EBED*60; %Packed bed contact time, sec 
SF = FLRT/(AREA*60); % Surface loading, cm/s
SFR = FLRT/(AREA*L); %Volumetric loading rate
%RE = 2*RHOL*RAD*SF/(VISC*EBED); 
RE = 2*RAD*RHOL*100*SF/(VISC*EBED); % Renold number
% End of fixed bed parameters

%% Calculate adsorption kinetics parameters
SPDFR = input(29,1);% Correction factor used in surface diffusivity calculation. Typical values range from 4 to 12

for i = 1:NCOMP
    % Calculation of liquid diffusivity based on Hayduk and Laudie
    % correlation, cm^2/s
    if DL(i)>=0
    else
        DL(i) = 13.26*10^(-5)/(VISC^(1.14)*VB(i)^(0.589));
    end
    
    % Calculation of Schmidt numbers
    SC(i) = VISC/(RHOL*DL(i)*100);
    
    % Calculation of film transfer coefficient based on Gnielinski
    % correlation, for fixed-bed, cm/s
    if KF(i) >= 0
    else
        KF(i) = (1+1.5*(1-EBED))*DL(i)*(2 + 0.644*RE^(1/2)*SC(i)^(1/3))/(2*RAD);
        % Calculate KF for preloaded GAC
        %KF(i)= 0.5*KF(i); % assuming preloaded GAC has reduced film transfer coefficient
    end
    
    % Calculation of surface diffusivity by equating the surface diffusion
    % flux to the pore diffusion flux with a correction factor
    if DS(i) >= 0
    else
        DS(i) = DL(i)*EPOR*CB0(i)*SPDFR/(TOR*RHOP*1000*XK(i)*CB0(i)^XN(i));
    end
    
    % Calculation of pore diffusivity from the liquid diffusivity
    if DP(i) >= 0
    else
        DP(i) = DL(i)/TOR;
    end  
end

CINF = 1;


%% Calculate dimensionless group for multiple compounds
QTE = 0.0;% Total equilibrium adsorbent phase conc., umole/g
DGT = 0.0;% Total solute distribution parameterD = DS./DP; 

for i = 1:NCOMP
    D(i) = DS(i)/DP(i); % D: ratio of surface to pore diffusivity, dimensionless
    QE(i) = XK(i)*CB0(i)^XN(i); % QE: adsorbent phase concentration in equilibrium with initial bulk phase concentration, umole/g 
    ST(i) = KF(i)*(1.0-EBED)*TAU/(EBED*RAD); % Modified Stanton number, dimensionless
    QTE = QTE + QE(i);
    DGS(i) = (RHOP*QE(i)*(1.0-EBED)*1000.0)/(EBED*CB0(i)); % DGS: surface solute distribution parameter, dimensionless
    DGP = EPOR*(1.0-EBED)/EBED; % DGP: pore solute distribution parameter, dimensionless
    DG(i)= DGS(i) + DGP; 
    DGT = DGT + DG(i); % Total solute distribution parameter, dimensionless
    EDS(i) = DS(i)*DGS(i)*TAU./RAD^2; % EDS: surface diffusion modulus, dimensionless
    EDP(i) = DP(i)*(DGP*TAU/RAD^2); % EDP: pore diffusion modulus, dimensionless
end

DGT1 = DGT+1;
STD = ST.*DGT1;
YM = QE./QTE; % Equilibrium adsorbent phase concentration fraction for each compound 
BIS = ST./EDS;
BIP = ST./EDP;
BVF = EBED*DGT; % Bed volumes fed to column

%% Convert independent variables and experimental data to dimensionless form
TCONV = 60/(TAU*(DGT+1)); % Independent variable (time) convertion constant, 1/min
TTOL = DTOL*TCONV; % total run time, dim
T0 = DT0*TCONV; % Initial time value in integration, dim

%% Orthogonal collocation
%Read in collocation constants
NC = input(2,1); % NC: number of radial collocation points, include 1 boundary = n+1
MC = input(3,1); % MC: number of axial collocation points, include 2 boundaries = n+2
ND = NC-1;
% MD = MC-1;
N1 = MC*(NC+1)-1;% N1: number of ordinary differential equations for each compounds
NEQ = N1*NCOMP;% NEQ: number of ordinary differential equations

[N1,ND,AZ,WR,BEDS,BEDP] = OrthogonalCollocation(NC,MC,ND,N1,NEQ,NCOMP,EDS,EDP,DGT1,DG,D);

%% Solve ordinary differential equations
% Use ode solver to integrate differential equations
tspan =  [T0,TTOL];%
Y0 = zeros(NEQ,1);% Initialize dependent variables
ODEdef_2 =@(t,Y) ODEdef(t,Y,NEQ,MC,NC,NCOMP,YM,Xn,N1,ND,CINF,AZ,DGT1,WR,BEDS,BEDP,STD,DG,XN); %Negative values in one column not good - same with zeroes
[t,Y] = ode15s(ODEdef_2, tspan, Y0);% Solve ODE equation sets using ODE solver

for i = 1:NCOMP
    breakthrough(:,i) = Y(:,N1*i); % effluent
end

time = t/TCONV; % operation time (min)
BV = time/EBCT;% bed volume passed through


%% Plot results

ind = find(breakthrough < 0, 1, 'last');
x = BV(ind+1:end);
y = breakthrough(ind+1:end);
aa = [x,y];

plot(x, y, 'linewidth', 1.5)
xlabel('Bed Volumes Treated')
ylabel('C/C_0')
set(gca, 'fontsize', 14)







