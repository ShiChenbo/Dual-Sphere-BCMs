f_begin=0.5;
f_end=7;
f_step=0.05;
len=round(1+(f_end-f_begin)/f_step);%频点数
signModes=1;%新模式的起始编号值
k=1;%初始化单独处理
f=f_begin+(k-1)*f_step
load([num2str(f*1e9),'.000000result.mat']);
NCMs=size(F,2);
ModesIndex=NaN./zeros(len,NCMs);%记录每个频点从大到小特征值的模式序号
for n=1:NCMs%电流单位化
    MS_n=abs(Val);
    
    if(MS_n(n)>=0)
        ModesIndex(k,n)= signModes;
        signModes=signModes+1;%更新新模式起始编号
    end
end
F_last=F;%只要关注的几个模式
Val_last=Val;

Qmatch=zeros(len,NCMs);%互相喜欢的模式
for k=2:len
    rmatch=logical(zeros(1,NCMs));%上一频点已经配对的模式
    cmatch=logical(zeros(1,NCMs));%这一频点已经配对的模式
    f=f_begin+(k-1)*f_step
    load([num2str(f*1e9),'.000000result.mat']);
    
    MS_n=abs((Val));
    for n=1:NCMs
        r1=abs(F_last'*F(:,n));
        [M,I]=max(r1);
        r0=abs(F_last(:,I)'*F);
        [M,Q]=max(r1);
        if(Q==n)%互相喜欢
            ModesIndex(k,n)=ModesIndex(k-1,I);
            rmatch(I)=1;
            cmatch(n)=1;
            F_last(:,I)=F_last(:,I)*NaN;
            %             Qmatch(k,n)=I;
        end
    end
    for n=1:NCMs
        if(cmatch(n)==0)
            cid=find(rmatch==0);%保存行（即上一频点）的序号
            if(~isempty(cid))
                rr=abs(F_last(:,cid)'*F(:,n));%返回列向量
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
    F_last=F;%只要关注的几个模式
    Val_last=Val;
end

% 根据跟踪结果转移每个频点每个模式的MS
MS=NaN./zeros(len,size(F,1));
Lam=NaN./zeros(len,size(F,1));
for k=1:len
    f=f_begin+(k-1)*f_step;
    load([num2str(f*1e9),'.000000result.mat']);
    for n=1:NCMs
        if(ModesIndex(k,n)>0)%存在这个模式
            MS(k,ModesIndex(k,n))=abs(Val(n));
            Lam(k,ModesIndex(k,n))=Val(n);
        end
    end
end

figure;
ka=2*pi/3e8*1e9*16.6667e-3*(f_begin:f_step:f_end);
plot(ka,MS(:,1:40),'linewidth',1.5)
xlim([0.41 2.3])
set(gca,'fontname','Times New Roman','fontsize',14)
grid on
xlabel('Frequency (GHz)')
ylabel('MS')
ylim([0 1])

figure;
plot((f_begin:f_step:f_end),MS(:,1:40),'linewidth',1.5)
set(gca,'fontname','Times New Roman','fontsize',14)
grid on
xlabel('Frequency (GHz)')
ylabel('MS')
ylim([0 1])