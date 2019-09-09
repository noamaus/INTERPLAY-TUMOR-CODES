%%%% Plots Figure 4B,D and E


load('DFIN11.mat') %% load gastrointestinal and endometrial data with MSI
load('SIG3.mat') %%  MSI-Aneuploidy signature

sm1=sum(DFIN11.MUT(SIG3,:));

MSI2=strcmp(DFIN11.MSI,'MSI-H');


canc=DFIN11.canc;
OS = DFIN11.OS;
death = DFIN11.death;
AN = DFIN11.AN;
MSI2=strcmp(DFIN11.MSI,'MSI-H');


subplot(2,3,1)
[pv2,h2] = ranksum(AN(sm1>median(sm1)),AN(sm1<=median(sm1)),'tail','left')
[pv3,h3] = ranksum(sm1(MSI2==1),sm1(MSI2==0),'tail','right');
[X,Y,T,AUC] = perfcurve(MSI2,sm1,1);
zz = sm1>0;


vv(zz==1) = {'High mutation rate'};
vv(zz==0) = {'Low mutation rate'};
[rr,pr] = corr(AN, sm1','type','Spearman')

%%% plot aneuploidy vs MSI-aneuploidy signature
scatter(AN, sm1,'filled')
hold on
legend(['rho = ',num2str(rr),'P = ',num2str(pr)])
xlabel('Aneuplody')
ylabel('signature mutation count')

subplot(2,3,2)
x1 = [OS(zz==1),death(zz==1)];
x2 = [OS(zz==0),death(zz==0)];
x1(isnan(x1))=0;
x2(isnan(x2))=0;

%%% plot KM for combined astrointestinal and endometrial tumors, based on
%%% MSI-aneuploidy signature load
[ppp] = logrank(x1,x2,'High mutation rate','Low mutation rate');
text(0.5,0.5,['P = ',num2str(ppp)])
text(4000,0.2,['n = ',num2str(length(x2))])
text(5000,0.8,['n = ',num2str(length(x1))])



%%% plot KM for each astrointestinal and endometrial tumor, based on
%%% MSI-aneuploidy signature load
uc=unique(DFIN11.canc);
for i = 1:length(uc)
    
    try
    subplot(2,3,i+2)
    canc1=strcmp(canc,uc(i));%%%
    OS1=OS(canc1);
    sm11=zz(canc1);
    death1=death(canc1);
    AN1=AN(canc1);
    MSI1=MSI2(canc1);
    zz1=zz(canc1)
     x1 = [OS1(zz1==1),death1(zz1==1)];
    x2 = [OS1(zz1==0),death1(zz1==0)];
    x1(isnan(x1))=0;
    x2(isnan(x2))=0;

    [p11(i)] = logrank(x1,x2,'High rate','Low rate');
    text(0.1,0.1,['P = ',num2str(p11(i))])
    title(uc(i))
    [pv11(i),h1] = ranksum(AN1(sm11>=median(sm11)),AN1(sm11<median(sm11)),'tail','left');
    [pv22(i),h2] = ranksum(MSI1(sm11>=median(sm11)),MSI1(sm11<median(sm11)),'tail','right');
    [X,Y,T,AUC1(i)] = perfcurve(MSI1,sm11,1);

    end
end