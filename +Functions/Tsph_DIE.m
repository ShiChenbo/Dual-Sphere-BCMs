function T = Tsph_DIE(varargin)
% This function is to calculate T matrix of a dielectric sphere
nargin_value = nargin; % number of inputs
switch nargin_value
    case 2
        deg = varargin{1};
        ka=varargin{2};
    case 3
        deg = varargin{1};
        freq = varargin{2};
        a= varargin{3};
        k=2*pi*freq*1e9/3e8;
        ka=k*a;
    case 4
        deg = varargin{1};
        freq = varargin{2};
        a= varargin{3};
        k=2*pi*freq*1e9/3e8;
        ka=k*a;
        epr=varargin{4};
        mu0=1;
        ep0=1;% assuming that the relative permittivity and permeability of free space are 1
        mur=1;% relative permeability of dielectric sphere
        gam_TE=mur/mu0;
        gam_TM=epr/ep0;
        kr=sqrt(gam_TE*gam_TM);% wave number in medium
    otherwise
        error('Invalid number of input arguments.');
end
%% Spherical Bessel Functions 
sph_Besselj=@(v,z)sqrt(pi./(2*z)).*besselj(v+0.5,z);
sph_Bessely=@(v,z)sqrt(pi./(2*z)).*bessely(v+0.5,z);
sph_Hankel=@(v,K,z)sqrt(pi./(2*z)).*besselh(v+0.5,K,z);

%% Differentiation of the Riccati-Bessel function
d_sph_R_j=@(v,z)(v+1).*sph_Besselj(v,z)-z.*sph_Besselj(v+1,z);
d_sph_R_y=@(v,z)(v+1).*sph_Bessely(v,z)-z.*sph_Bessely(v+1,z);
d_sph_R_h=@(v,K,z)(v+1).*sph_Hankel(v,K,z)-z.*sph_Hankel(v+1,K,z);

%% TE mode and TM mode. (tau=1 is TE mode, tau=2 is TM mode)
tTE=@(x,n)-(sph_Besselj(n,x).*d_sph_R_j(n,kr*x)-gam_TE*d_sph_R_j(n,x).*sph_Besselj(n,kr*x))...
    ./(sph_Hankel(n,2,x).*d_sph_R_j(n,kr*x)-gam_TE*d_sph_R_h(n,2,x).*sph_Besselj(n,kr*x));
tTM=@(x,n)-(sph_Besselj(n,x).*d_sph_R_j(n,kr*x)-gam_TM*d_sph_R_j(n,x).*sph_Besselj(n,kr*x))...
    ./(sph_Hankel(n,2,x).*d_sph_R_j(n,kr*x)-gam_TM*d_sph_R_h(n,2,x).*sph_Besselj(n,kr*x));
%% T Matrix
A=Functions.indexMatrix(deg);
N=size(A,2);
T=zeros(N);
for al=1:N
    l=A(1,al);
    if(A(4,al)==1)
        T(al,al)=tTE(ka,l);
    else
        T(al,al)=tTM(ka,l);
    end
end

