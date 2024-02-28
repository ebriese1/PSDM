function [N1,ND,AZ,WR,BEDS,BEDP] = OrthogonalCollocation(NC,MC,ND,N1,NEQ,NCOMP,EDS,EDP,DGT1,DG,D)

% spherical geometry,w = 1-w^2. Source: Finlayson 1972, Chapter 5 page 102
% in pdf    7x6 matrix
roots_radial = [0.654653670707977,1,0,0,0,0,0;
                0.468848793470714,0.830223896278567,1,0,0,0,0;
                0.363117463826178,0.677186279510738,0.899757995411460,1,0,0,0;...
                0.295758135586939,0.565235326996205,0.784483473663144,0.934001430408059,1,0,0;...
                0.249286930106240,0.482909821091336,0.686188469081758,0.846347564651872,...
                0.953309846642164,1,0;...
                0.215353955363794,0.420638054713673,0.606253205469846,0.763519689951815,...
                0.885082044222976,0.965245926503839,1];
%planar geometry, w = 1      9x9 matrix
roots_axial= [0.577350269100000,1,0,0,0,0,0,0,0;
              0.339981043500000,0.861136311500000,1,0,0,0,0,0,0;
              0.238619186000000,0.661209386400000,0.932469514200000,1,0,0,0,0,0;
              0.183434642400000,0.525532409900000,0.796666477400000,0.960289856400000,1,0,0,0,0;
              0.148874338900000,0.433395394100000,0.679409568200000,0.865063366600000,...
              0.973906528500000,1,0,0,0;
              0.125233408500000,0.367831498900000,0.587317954200000,0.769902674100000,...
              0.904117256300000,0.981560634200000,1,0,0;
              0.108054948700000,0.319112368900000,0.515248638300000,0.687292904800000,...
              0.827201315000000,0.928434883600000,0.986283808600000,1,0;
              0.0498749017000000,0.281603550700000,0.458016777600000,0.617876244400000,...
              0.755404408300000,0.865631202300000,0.944575023000000,0.989400934900000,1];

%Calculate Q, C, D and F matrices in Orthogonal Collocation based on roots
%of Legendre polynomials
Q_radial = zeros(NC,NC); Q_axial = zeros(MC,MC);
C_radial = zeros(NC,NC); C_axial = zeros(MC,MC);
D_radial = zeros(NC,NC); D_axial = zeros(MC,MC);
F_radial = [1:NC]; F_axial = [1:MC];

% a = 3, for spherical geometry
for j = 1:NC
   for i = 1:NC
       Q_radial (j,i) =  roots_radial(NC-1,j)^(2*i-2);
       C_radial (j,i) =  (2*i-2)*roots_radial(NC-1,j)^(2*i-3);
       D_radial (j,i) =  (2*i-2)*(2*i-3)*roots_radial(NC-1,j)^(2*i-4);
       F_radial (i) = 1/(2*i+1);
   end
end

% a=1, for planar geometry
for j = 1:MC
   for i = 1:MC
       Q_axial (j,i) =  roots_axial(MC-1,j)^(2*i-2);
       C_axial (j,i) =  (2*i-2)*roots_axial(MC-1,j)^(2*i-3);
       D_axial (j,i) =  (2*i-2)*(2*i-3)*roots_axial(MC-1,j)^(2*i-4);
       F_axial (i) = 1/(2*i-1);
   end
end

%Calculate A, B and W matrices in Orthogonal Collocation 
BR = D_radial/Q_radial;
WR = F_radial/Q_radial;
AZ = C_axial/Q_axial;
BEDS = zeros (NCOMP,ND,NC);
BEDP = zeros (NCOMP,ND,NC);

for i = 1:NCOMP    
    for j =1:ND
        for k = 1:NC
            BEDS(i,j,k) = (EDS(i)+D(i)*EDP(i))*DGT1/DG(i)*BR(j,k);
            BEDP(i,j,k) = EDP(i)*(1-D(i))*DGT1/DG(i)*BR(j,k);
        end
    end
end 

end

