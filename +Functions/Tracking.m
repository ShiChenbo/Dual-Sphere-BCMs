function Tracking(DataFileName,fbeg,fend,fstp,radius)
len=round(1+(fend-fbeg)/fstp);% number of frequencies
signModes=1;% Mode index
k=1;
f=fbeg;
load([DataFileName '\' num2str(f*1e9),'.000000result.mat']);
NCMs=size(F,2);
ModesIndex=NaN./zeros(len,NCMs);% Tracking info
%% initialize
for n=1:NCMs
    MS_n=abs(Val);    
    if(MS_n(n)>=0)
        ModesIndex(k,n)= signModes;
        signModes=signModes+1;
    end
end
F_last=F;
Val_last=Val;
%% Begin Tracking
for k=2:len
    rmatch=logical(zeros(1,NCMs));
    cmatch=logical(zeros(1,NCMs));
    f=fbeg+(k-1)*fstp;
    load([DataFileName '\' num2str(f*1e9),'.000000result.mat']);    
    MS_n=abs((Val));
    for n=1:NCMs
        r1=abs(F_last'*F(:,n));% first group of correlation coefficient
        [M,I]=max(r1);
        r0=abs(F_last(:,I)'*F);% second group of correlation coefficient
        [M,Q]=max(r1);
        if(Q==n)% means that the max correlation match each group
            ModesIndex(k,n)=ModesIndex(k-1,I);
            rmatch(I)=1;
            cmatch(n)=1;
            F_last(:,I)=F_last(:,I)*NaN;           
        end
    end
    for n=1:NCMs % the rest untracked modes using max correlation coefficient
        if(cmatch(n)==0)
            cid=find(rmatch==0);
            if(~isempty(cid))
                rr=abs(F_last(:,cid)'*F(:,n));
                [M,I]=max(rr);
                if(M>=0)
                    ModesIndex(k,n)=ModesIndex(k-1,cid(I));
                    rmatch(cid(I))=1;
                    cmatch(n)=1;
                    F_last(:,cid(I))=F_last(:,cid(I))*NaN;
                end
            end
        end
    end
    F_last=F;
    Val_last=Val;
end

%% Sort the tracked data
MS=NaN./zeros(len,size(F,1));
Lam=NaN./zeros(len,size(F,1));
for k=1:len
    f=fbeg+(k-1)*fstp;
    load([DataFileName '\' num2str(f*1e9),'.000000result.mat']);
    for n=1:NCMs
        if(ModesIndex(k,n)>0)
            MS(k,ModesIndex(k,n))=abs(Val(n));
            Lam(k,ModesIndex(k,n))=Val(n);
        end
    end
end
%% plot the results
figure;
ka=2*pi/3e8*1e9*radius*(fbeg:fstp:fend);
plot(ka,MS,'linewidth',1.5)
xlim([min(ka) max(ka)])
set(gca,'fontname','Times New Roman','fontsize',8)
grid on
xlabel('$k\Sigma$','interp','latex')
ylabel('MS')
ylim([0 1])

figure;
plot((fbeg:fstp:fend),MS,'linewidth',1.5)
xlim([fbeg fend])
set(gca,'fontname','Times New Roman','fontsize',8)
grid on
xlabel('Frequency (GHz)')
ylabel('MS')
ylim([0 1])
end

