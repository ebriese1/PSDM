%% Function: ODE equation sets of PSDM (including boundary conditions)
function [dYdt] = MCS_PSDM_ODEdef(t,Y, NEQ,MC,NC,NCOMP,YM,Xn,N1,ND,CINF,AZ,DGT1,WR,BEDS,BEDP,STD,DG,XN)
% Initialize dependent variables
dYdt = zeros(NEQ,1);
Cp = zeros(NEQ,1); %initialize liquid phase concentration in pores

% Calculate Cp at each radial and axial position using IAST
% Initialize constants
ii = 0;
jj = 0;
% Calculate Cp
for j = 1: MC
    for k = 1: NC
        QTE = 0; 
        YT0 = 0;
            for i = 1:NCOMP
                ii = ii +1;
                Z(i) = YM(i)*Y(ii);
                QTE = QTE + Z (i);
                YT0 = YT0 + Xn(i)*Z(i);
                ii = ii + N1-1;
            end
        
        for i = 1:NCOMP
            jj = jj+1;
            if ((QTE <= 0.0) || (YT0 <= 0.0))
                Cp(jj) = 0;
            else
                Z(i) = Z(i)./QTE;
                Q0(i) = YT0 .* XN(i)./YM(i);
                if ( Xn(i) * log10(Q0(i)) < -20.0)
                    Cp(jj) = 0.0;
                else
                    Cp(jj) = Z(i)*Q0(i)^Xn(i);
                end
            end
            jj = jj+N1-1;
        end
            if k < ND
                ii = (j-1)*ND+k;
                jj = (j-1)*ND+k;
            else
                ii = (j-1)+MC*ND;
                jj = (j-1)+MC*ND;
            end
    end
    ii = ND * j;
    jj = ND * j;
end

%%%%%% portion above is good


for i = 1: NCOMP
    for k = 2: MC
        if (Cp((i-1)*N1+MC*ND+k) <= 0)
            CBS(i,k) = STD(i)*Y((i-1)*N1+MC*(NC-1)+(MC-1)+k);
        else
            CBS(i,k) = STD(i)*(Y((i-1)*N1+MC*(NC-1)+(MC-1)+k)-Cp((i-1)*N1+MC*(NC-1)+k));
        end
    end
    
        WrY = zeros(NCOMP,MC);
        AzC = zeros(NCOMP,MC);
        BrY = zeros(NCOMP,ND,MC);
        BrCp = zeros(NCOMP,ND,MC);
        
    for k = 1:MC   
        for j =  1: ND
            for m = 1: ND
                BrY(i,j,k) = BEDS(i,j,m) * Y(((i-1)*N1+(k-1)*(NC-1))+m) + BrY(i,j,k);
                BrCp(i,j,k) = BEDP(i,j,m) * Cp(((i-1)*N1+(k-1)*(NC-1))+m) + BrCp(i,j,k);
            end
            BrY(i,j,k) = BEDS(i,j,NC) * Y((i-1)*N1+MC*(NC-1)+k)+BrY(i,j,k);
            BrCp(i,j,k) = BEDP(i,j,NC)* Cp((i-1)*N1+MC*(NC-1)+k)+BrCp(i,j,k);
        end
        
     % Intraparticle phase mass balance (excluding boundary)
        for j = 1: ND   
            dYdt (((i-1)*N1+(k-1)*(NC-1)+j)) = BrY(i,j,k)+BrCp(i,j,k);
            WrY(i,k) = WR(j)*dYdt(((i-1)*N1+(k-1)*(NC-1)+j))+WrY(i,k);
        end
    end
    
        % Boundary conditions
        % Liquid-Solid boundary layer mass balance at column entrance
          dYdt((i-1)*N1+MC*(NC-1)+1) = (STD(i)/DG(i)*(CINF-Cp((i-1)*N1+MC*(NC-1)+1)) - WrY(i,1))/WR(NC); % WrY(i,1) 
          
    for k = 2: MC
           % Liquid-Solid boundary layer mass balance within column
           dYdt((i-1)*N1+MC*ND+k) = (CBS(i,k)/DG(i)- WrY(i,k))/WR(NC);
           
           for p = 2: MC
               AzC(i,k) = AZ(k,p)*Y((i-1)*N1+MC*(NC-1)+(MC-1)+p) + AzC(i,k);
           end
           
            % Liquid phase mass balance
           dYdt ((i-1)*N1+MC*(NC-1)+(MC-1)+k) = DGT1*(-AZ(k,1)*CINF-AzC(i,k))-3*CBS(i,k);
           
    end 
       
end
    
end