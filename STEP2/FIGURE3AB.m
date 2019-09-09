%% plot Figure 3A-B

load('DFIN1.mat') %%gastrointestinal and endometrial data
load('DFIN2.mat') %%all other tumors data
load('SIG1.mat') %%Repair signature
load('SIG2.mat') %%Apoptosis signature
load('DFIN.mat') %%PAN CANCER data

 
%% Get signatures load
sc1 = (sum(DFIN.MUT(SIG1,:)));%%Repair
sc2 = (sum(DFIN.MUT(SIG2,:)));%%APOP
scc = (sum(DFIN.MUT(SIG2,:))+1)./(sum(DFIN.MUT(SIG1,:))+1);%%APOP/repair


load('UC1.mat')
RMR=ones(1,length(UC1));
PMR=ones(1,length(UC1));
RMRR=ones(1,length(UC1));
RMRA=ones(1,length(UC1));
PMRR=ones(1,length(UC1));
PMRA=ones(1,length(UC1));


P1=[];ANL=[];ANH=[];LL=[];LH=[];
APOP1=[];REP1=[];APL1=[];RPL1=[];
[RM,PM] = corr(DFIN.AN,scc','type','Spearman');


for i = 1:length(UC1)
    i
%     subplot(3,16,i)
    cc = find(strcmp(DFIN.canc,UC1(i)));
    s2 = scc(cc);
    sr = sc2(cc);%%APOP
    sa = sc1(cc);
    AN1 = DFIN.AN(cc);

    
    [RMR(i),PMR(i)] = corr(AN1,s2','type','Spearman'); %%% Correlate ratio

    [RMRR(i),PMRR(i)] = corr(AN1,sr','type','Spearman');%%Correlate apoptosis

    [RMRA(i),PMRA(i)] = corr(AN1,sa','type','Spearman');%%Correlate repair

    [P1(i),H1(i)] = ranksum(AN1(s2<=1),AN1(s2>1)','tail','left');
    
    
    %%%DEL
    CNT1(i) = nnz(s2<1);%%Higure in REPAIR
    CNT2(i) = nnz(s2>1);%%Higure in APOP
    APOP1 = [APOP1;sr'];
    APL1 = [APL1;i*ones(size(sr'))];
    REP1 = [REP1;sa'];
    RPL1 = [RPL1;i*ones(size(sa'))];
    
    
    ANL = [ANL;AN1(s2<1)];
    ANH = [ANH;AN1(s2>1)];
    LL = [LL;i*ones(size(AN1(s2<1)))];
    LH = [LH;i*ones(size(AN1(s2>1)))];
end

subplot(3,1,1)
boxplot(ANL,LL,'Symbol','g.','Jitter',0.1)

ylim([-1,40])

subplot(3,1,2)
boxplot(ANH,LH,'Symbol','g.','Jitter',0.1)
ylim([-1,40])



subplot(3,1,3)
xx=[RMR;RMRR;RMRA]
bar(xx',1.3)
legend('Ratio','Apoptosis','Repair')


xp1=find(PMR<=0.1);
yp1=RMR(xp1);
yp1(yp1>=0) = yp1(yp1>=0)+0.025;
yp1(yp1<0) = yp1(yp1<0)-0.125;
xp1=xp1-0.35;
tt(xp1>0) = {'*'};
text(xp1,yp1,tt,'Color','blue','FontSize',25)

xp2=find(PMRR<=0.1);
yp2=RMRR(xp2);
yp2(yp2>=0) = yp2(yp2>=0)+0.025;
yp2(yp2<0) = yp2(yp2<0)-0.125;
tt2(xp2>0) = {'*'};
text(xp2,yp2,tt2,'Color','red','FontSize',25)

xp3=find(PMRA<=0.1);
yp3=RMRA(xp3);
yp3(yp3>=0) = yp3(yp3>=0)+0.025;
yp3(yp3<0) = yp3(yp3<0)-0.125;
tt3(xp3>0) = {'*'};
text(xp3,yp3,tt3,'Color','yellow','FontSize',25)


