function T = Tsph_Layer(varargin)
% This function is used to compute the T-matrix of a multilayer sphere. The case in the paper is just a special case where the kernel is PEC
nargin_value = nargin; % number of inputs
switch nargin_value
    case 2
        deg = varargin{1};
        ka=varargin{2};
    case 3
        deg = varargin{1};
        freq = varargin{2};
        r= varargin{3};% radii of different layers
        k=2*pi*freq*1e9/3e8;
        ka=k*r;
    case 4
        deg = varargin{1};
        freq = varargin{2};
        r= varargin{3};
        k0=2*pi*freq*1e9/3e8;
        er=varargin{4};% relative permittivity of different layers
        er(end+1)=1;
        r(end+1)=inf;
        Nlay=length(r);% number of layers
        mr=ones(Nlay,1);% relative permeability of different layers
    otherwise
        error('Invalid number of input arguments.');
end
%% Spherical Bessel Functions 
sph_Besselj=@(v,z)sqrt(pi./(2*z)).*besselj(v+0.5,z);
sph_Bessely=@(v,z)sqrt(pi./(2*z)).*bessely(v+0.5,z);
sph_Hankel=@(v,K,z)sqrt(pi./(2*z)).*besselh(v+0.5,K,z);

%% Riccati-Bessel Functions
R_j=@(v,z)z.*sph_Besselj(v,z);
R_h=@(v,K,z)z.*sph_Hankel(v,K,z);
%% Differentiation of the Riccati-Bessel function
dR_j=@(v,z)(v+1).*sph_Besselj(v,z)-z.*sph_Besselj(v+1,z);
dR_h=@(v,K,z)(v+1).*sph_Hankel(v,K,z)-z.*sph_Hankel(v+1,K,z);
%% T Matrix

for ilay=1:Nlay-1
    gamma1(ilay)=sqrt(mr(ilay)/er(ilay))/sqrt(mr(ilay+1)/er(ilay+1));
    gamma2(ilay)=sqrt(mr(ilay+1)/er(ilay+1))/sqrt(mr(ilay)/er(ilay));
end

A=Functions.indexMatrix(deg);
N=size(A,2);
T=zeros(N);
for al=1:N
    l=A(1,al);
    k0=2*pi*freq*1e9/3e8;
    for ilay=1:Nlay
        k(ilay)=sqrt(er(ilay)*mr(ilay))*k0;
    end
    if(A(4,al)==1)
        if(isinf(er(1)))% Determine whether the inner layer is PEC
            % PEC case recursive start from 2
            t(1)=-R_j(l,k(2)*r(1))/R_h(l,2,k(2)*r(1));
            for ilay=2:Nlay-1
                a(ilay)=gamma1(ilay)*R_h(l,2,k(ilay)*r(ilay))*dR_j(l,k(ilay+1)*r(ilay))-dR_h(l,2,k(ilay)*r(ilay))*R_j(l,k(ilay+1)*r(ilay));
                b(ilay)=gamma1(ilay)*R_j(l,k(ilay)*r(ilay))*dR_j(l,k(ilay+1)*r(ilay))-dR_j(l,k(ilay)*r(ilay))*R_j(l,k(ilay+1)*r(ilay));
                c(ilay)=gamma1(ilay)*R_h(l,2,k(ilay)*r(ilay))*dR_h(l,2,k(ilay+1)*r(ilay))-dR_h(l,2,k(ilay)*r(ilay))*R_h(l,2,k(ilay+1)*r(ilay));
                d(ilay)=gamma1(ilay)*R_j(l,k(ilay)*r(ilay))*dR_h(l,2,k(ilay+1)*r(ilay))-dR_j(l,k(ilay)*r(ilay))*R_h(l,2,k(ilay+1)*r(ilay));
                t(ilay)=-(a(ilay)*t(ilay-1)+b(ilay))/(c(ilay)*t(ilay-1)+d(ilay));
            end
            T(al,al)=t(Nlay-1);
        end
        
    else
        if(isinf(er(1)))% Determine whether the inner layer is PEC
            t(1)=-dR_j(l,k(2)*r(1))/dR_h(l,2,k(2)*r(1));
            for ilay=2:Nlay-1
                a(ilay)=gamma2(ilay)*R_h(l,2,k(ilay)*r(ilay))*dR_j(l,k(ilay+1)*r(ilay))-dR_h(l,2,k(ilay)*r(ilay))*R_j(l,k(ilay+1)*r(ilay));
                b(ilay)=gamma2(ilay)*R_j(l,k(ilay)*r(ilay))*dR_j(l,k(ilay+1)*r(ilay))-dR_j(l,k(ilay)*r(ilay))*R_j(l,k(ilay+1)*r(ilay));
                c(ilay)=gamma2(ilay)*R_h(l,2,k(ilay)*r(ilay))*dR_h(l,2,k(ilay+1)*r(ilay))-dR_h(l,2,k(ilay)*r(ilay))*R_h(l,2,k(ilay+1)*r(ilay));
                d(ilay)=gamma2(ilay)*R_j(l,k(ilay)*r(ilay))*dR_h(l,2,k(ilay+1)*r(ilay))-dR_j(l,k(ilay)*r(ilay))*R_h(l,2,k(ilay+1)*r(ilay));
                t(ilay)=-(a(ilay)*t(ilay-1)+b(ilay))/(c(ilay)*t(ilay-1)+d(ilay));
            end
            T(al,al)=t(Nlay-1);
        end
    end
end
end

