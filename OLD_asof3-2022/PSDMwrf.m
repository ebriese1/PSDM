%% Read in data 
clear;clc;

% Read data from excel spreadsheet "PSDM.xlsx"
input = xlsread('PSDMwrf.xlsx');
inputRawData = xlsread('PSDMwrf.xlsx', "RawData");
RawDataBVT1 = inputRawData(1:end,1); % pH = 8.3 Test 1
RawDataC1 = inputRawData(1:end,2);

% RawDataBVT2 = inputRawData(1:end,3); % pH = 6.0 Test 2 
% RawDataC2 = inputRawData(1:end,4);
% 
% RawDataBVT3 = inputRawData(1:end,5); % pH = 8.3 
% RawDataC3 = inputRawData(1:end,6);
% 
% RawDataBVT4 = inputRawData(1:end,7); % pH = 8.3 
% RawDataC4 = inputRawData(1:end,8);

% Read index No. for multiruns
% TestID = input(4,3);

% Read data for model input
% input(1,1)=2;
% input(6,2)=10;
% input(7,2)=0.4;
% input(8,2)=0.1;
% input(9,2)=50; % doesnt matter
% input(10,2)=input(10,1);

NCOMP = input(1,1); % NCOMP: number of components
NCOMP = 1;
trial=4;
XK = input(6,1:NCOMP); %Freundlich isotherm capacity K,(mg/g)(L/mg)^(1/n)
%XK = 1.98; % 1.98 pH = 8.3 and 2.64 pH = 6.0         K,(ug/mg)(L/ug)^(1/n)
%XK = 2.64;
%XK = 14; % pH=8.3 K,(mg/g)(L/mg)^(1/n)
XK = 24; % pH=6.0


%XN = input(7,1:NCOMP); % Freundlich isotherm exponents 1/n, dimensionless
%XN = 0.36; % 0.36 pH = 8.3 and 0.33 pH = 6.0
XN = 0.3;
Xn = 1./XN; %Freundlich isotherm exponent reciprocal n, dim.
CB0 = input(8,1:NCOMP);% Initial concentrations, mg/L
% CB0 = input(8,1:NCOMP)*1000;% Initial concentrations, microg/L
MW = input(9,1:NCOMP);% Molecular weight, g/mol
VB = input(10,1:NCOMP);% molar volume at the boiling point
VB = 30.5; % For As in cm^3/mol

L = input(13,1);% L: length of column,cm
DIA = input(14,1); % diameter of column, cm
WT = input(15,1); % WT: mass of adsorbent in column, grams
FLRT = input(16,1); % Flow rate, ml/min
TEMP = input(18,1); % water temperature, celsius degree; 
TEMP = 25;
TEMP = TEMP + 273.15; % water temperature, Kelvin;
TEMPA = TEMP/324.65; % Absolute water temperature, dim.;

RHOP = input(21,1);% apparent adsorbent density, g/cm^3
EPOR = input(22,1);% Particle porosity: fraction of volumetric space in adsorbent phase unoccupied by adsorbent (pore void fraction)
RAD = input(23,1); % Radius of adsorbent, cm
% TOR = input(24,1); % tortuosity of particle
TOR = (2-EPOR).^2./EPOR;

% Density of water at a certain temperature, g/cm^3 = 1000 kg/m^3
% Calculated using a correlation developed at Michigan Technological University
RHOL = (-1.41768+8.976651.*TEMPA-12.2755.*TEMPA.^2+7.45844.*TEMPA.^3-1.738492.*TEMPA.^4)*0.98396; % [g/cm^3]

% Viscosity of water at a certain temperature,cP (=0.001 kg/m/s)
% Calculated using a correlation shown by Reid (1987) and Yaws (1976)
VISC = exp(-24.71+4209./TEMP+0.04527.*TEMP-3.376.*10.^(-5).*TEMP.^2); % 1 cP; where 1000 cP = 1[kg/(m*s)]


KF = input(26,1:NCOMP); % Film transfer coefficient, cm/s
DS = input(27,1:NCOMP); % Surface diffusivity,cm^2/s
DP = input(28,1:NCOMP);% Pore diffusivity based on pore void fraction,cm^2/s
DL = input(34,1:NCOMP); % Liquid diffusivity, cm^2/s
SC = input(35,1:NCOMP); % Schmidt number, dim.

Days=10.5; %10.5 AND 5.5
DTOL = input(31,1); % total run time, min.
DTOL = Days*24*60;
DTOL = 10000;
DSTEP = 1; % Time step added in integration, min.
DT0 = 0.1; % Initial time value in integration, min.

%% Calculate fixed bed parameters
AREA = pi*DIA*DIA/4; % Cross sectional area of differential slice in bed, cm^2
BEDVOL = 1000; % Bed volume, ml
% WT = EPOR.*RHOP.*pi().*DIA.^2/4.*L;
% WT = WT/2;
EBED = 1.0 - WT/(BEDVOL*RHOP); % Bed porosity: void fraction (Fraction of volumetric space in reactor unoccupied by adsorbent)

% EBED = 0.428;
EBCT = BEDVOL/FLRT/EBED;  % Empty bed contact time, min
TAU = EBCT*EBED*60; %Packed bed contact time, sec 
SF = FLRT/(AREA*60); % Surface loading, cm/s
InterStitVelo = FLRT/AREA/EBED/60; % Interstitial Velocity cm/s
%RE = 2*RHOL*RAD*SF/(VISC*EBED); 
phi=1.25; % Particle Shape Factor (from Kiril
% RE = phi*2*RAD*RHOL*SF/(VISC*EBED); % Reynold's number, /100^2 for converting to m
% RE = phi*2*RAD*RHOL*SF/(VISC*EBED)/(100^2); % Reynold's number, /100^2 for converting to m
RE = phi*2*RAD*(RHOL)*InterStitVelo/(VISC*0.1); % this is correct
% RE = DIA*(RHOL)*InterStitVelo/(VISC*0.01); % this is correct
% End of fixed bed parameters


% DTOL = 100000/(FLRT/BEDVOL);



%% Calculate adsorption kinetics parameters
% somenumber=8
% [changeavariable]=linspace



% for iterations=1:somenumber

SPDFR = input(29,1);% Correction factor used in surface diffusivity calculation. Typical values range from 4 to 12
SPDFR = 8;

DL = nan;
% SC = nan; % Always calc'd unless given - usually not given
KF = nan;
DS = nan;
% DP = nan;

VB = 113;

for i = 1:NCOMP
    % Calculation of liquid diffusivity based on Hayduk and Laudie
    % correlation, cm^2/s
    if DL(i)>=0
    else
        DL(i) = 13.26*10^(-5)/(VISC^(1.14)*VB(i)^(0.589)); % [cm^2/s]
    end
    
    % Calculation of Schmidt numbers
    SC(i) = (VISC*0.01)/(RHOL*DL(i));
    
    % Calculation of film transfer coefficient based on Gnielinski
    % correlation, for fixed-bed, cm/s
    if KF(i) >= 0
    else
        KF(i) = ((1+1.5*(1-EBED)) * DL(i)) / (2*RAD) * (2 + 0.644*RE^(1/2)*SC(i)^(1/3)) ;  % [cm/s]
        % Calculate KF for preloaded GAC
        %KF(i)= 0.5*KF(i); % assuming preloaded GAC has reduced film transfer coefficient
    end
    
    % Calculation of surface diffusivity by equating the surface diffusion
    % flux to the pore diffusion flux with a correction factor
    if DS(i) >= 0
    else
         DS(i) = DL(i)*EPOR*(CB0(i)./1000)*SPDFR / (TOR*RHOP*XK(i)*CB0(i)^XN(i)); % when mg/L used
%        DS(i) = DL(i)*EPOR*(CB0(i)./1000)*SPDFR / (TOR*RHOP*1000*XK(i)*CB0(i)^XN(i)); % when microg/L used
    end
    
    % Calculation of pore diffusivity from the liquid diffusivity
    if DP(i) >= 0
    else
        DP(i) = DL(i)/TOR;
    end  
end

qe = XK*CB0^(XN);
Dg = RHOP*qe*(1-EBED)/(EBED*CB0/1000);
Biot = KF*RAD*(1-EBED)/(DS*Dg*EBED);
Stanton = KF*TAU*(1-EBED)/(EBED*RAD);

        KF2 = 1.85 * ((1-EBED)/EBED)^(1/3) * RE^(1/3) * SC^(1/3);
zz= [DL;SC;KF;DS;DP];

MainParametersNAMES= {'XK'; 'XN'; 'KF'; 'DS'; 'RAD'; 'FLRT'; 'L'; 'DIA'; 'WT'; 'TEMP'; 'EBED'; 'RHOP'; 'EPOR'; 'phi'; 'TOR'};
MainParametersVALUES= [XK; XN; KF; DS; RAD; FLRT; L; DIA; WT; TEMP; EBED; RHOP; EPOR; phi; TOR];
MainParametersTable = array2table(MainParametersVALUES,'RowNames',MainParametersNAMES);

CINF = 1; 

%% Up to this point changing XK, XN, or CB0 will not change matrix sizes


%% Calculate dimensionless group for multiple compounds
QTE = 0.0;% Total equilibrium adsorbent phase conc., mg/g
DGT = 0.0;% Total solute distribution parameterD = DS./DP; 

for i = 1:NCOMP
    D(i) = DS(i)./DP(i); % D: ratio of surface to pore diffusivity, dimensionless
    QE(i) = XK(i).*CB0(i).^XN(i); % QE: adsorbent phase concentration in equilibrium with initial bulk phase concentration, mg/g 
    ST(i) = KF(i).*(1.0-EBED).*TAU./(EBED.*RAD); % Modified Stanton number, dimensionless
    QTE = QTE + QE(i);
         DGS(i) = (RHOP.*QE(i).*(1.0-EBED).*1000.0)./(EBED.*CB0(i)); % DGS: surface solute distribution parameter, dimensionless
%        DGS(i) = (RHOP.*1000.*QE(i).*(1.0-EBED).*1000.0)./(EBED.*CB0(i)); % when microg/L used
    DGP = EPOR.*(1.0-EBED)./EBED; % DGP: pore solute distribution parameter, dimensionless
    DG(i)= DGS(i) + DGP; 
    DGT = DGT + DG(i); % Total solute distribution parameter, dimensionless
    EDS(i) = DS(i).*DGS(i).*TAU./RAD.^2; % EDS: surface diffusion modulus, dimensionless
    EDP(i) = DP(i).*DGP(1).*TAU./RAD.^2; % EDP: pore diffusion modulus, dimensionless
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

NC = 5;
MC = 9;
ND = NC-1;
% MD = MC-1;
N1 = MC*(NC+1)-1;% N1: number of ordinary differential equations for each compounds
NEQ = N1*NCOMP;% NEQ: number of ordinary differential equations

[N1,ND,AZ,WR,BEDS,BEDP] = OrthogonalCollocation(NC,MC,ND,N1,NEQ,NCOMP,EDS,EDP,DGT1,DG,D);

%% Solve ordinary differential equations
% Use ode solver to integrate differential equations
tspan =  [T0,TTOL];%
Y0 = zeros(NEQ,1);% Initialize dependent variables
ODEdef_2 =@(t,Y) ODEdef(t,Y,NEQ,MC,NC,NCOMP,YM,Xn,N1,ND,CINF,AZ,DGT1,WR,BEDS,BEDP,STD,DG,XN);
[t,Y] = ode15s(ODEdef_2, tspan, Y0);% Solve ODE equation sets using ODE solver

for i = 1:NCOMP
    breakthrough(:,i) = Y(:,N1*i); % effluent
end

time = t/TCONV; % operation time (min)
BV = time/EBCT;% bed volume passed through


%% Plot results
% STOREx=nan(1000,NCOMP);
% STOREy=nan(1000,NCOMP);

% for i=1:NCOMP
ind = find(breakthrough(:,i) < 0, 1, 'last');
x = BV(ind+1:end,1);
y = breakthrough(ind+1:end,i);

a=[x y];

% % % CASE START
% breakthrough(1:ind,1)=0;
% x = time./EBCT;
% y = breakthrough;
% % CASE END

% % % % CASE START
% breakthrough(1:ind,1)=0;
% xxx=time(ind,1);
% xx=time(ind+1:end);
% x=xx-xxx;
% x = x./EBCT;
% y = breakthrough(ind+1:end);
% % % CASE END


% time(1:ind,1) = 0;
% breakthrough(1:ind,1) = 0;
% x=time(ind:end);
% y=breakthrough(ind:end);

% x = x./60./24;
% x=BV;
% BV = x*FLRT/(BEDVOL*EPOR);
% ttime = x.*TCONV;
figure() 
plot(x,y, 'k', 'linewidth', 3, 'DisplayName','Model Curve');

figure()
hold on
% plot(ttime,y,'linewidth',2)
p1=plot(x,y, 'k', 'linewidth', 3, 'DisplayName','Model Curve');
% p1=plot(x,y, 'linewidth', 1, 'DisplayName','Experimental Parameters');
% p1=plot(x,y, 'linewidth', 1, 'DisplayName','Calculating Ds Only');

% % p2=plot(RawDataBVT2, RawDataC2,'or','linewidth', 2, 'DisplayName','Raw Data (pH=8.3)');
% % % p3=plot(RawDataBVT2, RawDataC2,'or','linewidth', 2, 'DisplayName','Raw Data (pH=6.0)');
% % % title('WRF Raw Data vs. Model Prediction (Arsenic and E33)')
% % title('As Removal by E33 at pH=8.3','FontSize',28, 'Color', 'r', 'FontWeight', 'bold');
% % % xlabel('Time (day)')
% % xlabel('BVT')
% % ylabel('Removal Efficiency (%)')
% % ax = gca;
% % ax.XAxis.FontSize=24;
% % ax.YAxis.FontSize=24;
% % 
% % lgd=legend([p1,p2],'Location','SouthEast');
% % lgd.FontSize=14;
% % hold off
% % 
% % FinalParametersNAMES= {'min BVT'; 'max BVT'; 'min Conc'; 'max Conc'; 'min exp BVT'; 'max exp BVT'; 'min exp Conc'; 'max exp Conc'};
% % FinalParametersVALUES= [min(x); max(x); min(y); max(y); min(RawDataBVT2); max(RawDataBVT2); min(RawDataC2); max(RawDataC2); ];
% % FinalParametersTable = array2table(FinalParametersVALUES,'RowNames',FinalParametersNAMES);
% % 
% % % end %the end of changing 
% % 
% % % ind = find(breakthrough(:,2) < 0, 1, 'last');
% % % x = BV(ind+1:end);
% % % y = breakthrough(ind+1:end,2);
% % 
% % % % % ---- Good ---- % % % 
% % % for i = 1:NCOMP
% % %     ind(1,i) = find(breakthrough(:,i) < 0, 1, 'last');
% % % end
% % % for i = 1:NCOMP
% % % x{i}= BV(ind(i)+1:end);
% % % y{i} = breakthrough(ind(i)+1:end,i);
% % % STORE{i}=[x{i},y{i}];
% % % end
% % % 
% % % 
% % % figure()
% % % hold on
% % % BiSol=cellfun(@plot,x,y);
% % % BiSol(1).Marker = 'o';
% % % % BiSol(2).Marker = '+';
% % % % % ---- Good ---- % % % hold off
% % 
% % ZZZZ=[x y];
% % 
% % % % 
% % % plot(STOREx,STOREy, 'linewidth', 1.5)
% % % xlabel('Bed Volumes Treated')
% % % ylabel('C/C_0')
% % % set(gca, 'fontsize', 14)
% % % ylim([0 1])
% % trialdj=cell2mat(STORE);
% % % arrayfun(@(x) plot(STORE(1,1:NCOMP),STORE(2,1:NCOMP),STORE))
% % cellplot(STORE(3,:),STORE(4,:))
% % arrayfun(@(x) plot(x{:}(:,1),x{:}(:,2)),data)
% % 
% % plot(x(:,:),y(:,:), 'linewidth', 1.5)
% % xlabel('Bed Volumes Treated')
% % ylabel('C/C0 (%)')
% % set(gca, 'fontsize', 14)
% % % ylim([0 1])
% % 
% % 
% % % STOREx(1:length(x),iter)=x;
% % % STOREy(1:length(x),iter)=y;
% % 
% % 
% % 
% % 
