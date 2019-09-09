%%%Plot Figure 1B

%%% Load data parted to gastrointestinal and endo/All other
load('DFIN1.mat') %%gastrointestinal and endometrian data
load('DFIN2.mat') %%all other tumors data
load('PANDRIVE.mat') %%Pan cancer drivers
cnt=1;
ps1=[];ps2=[];dAUC1=[];dAUC2=[];

[a,b,c] = intersect(DFIN.gene,PANDRIVE);


%%%Plot KM analysis for each thershold of #drivers, for each tumor cluster

for i = 1:2:8

THR=i;

canc=DFIN2.canc;
OS = DFIN2.OS;
death = DFIN2.death;
DRC = sum(DFIN2.MUT(b,:));
ss1 = DRC>THR;


subplot(2,4,cnt)
x1 = [OS(ss1==1),death(ss1==1)];
x2 = [OS(ss1==0),death(ss1==0)];
x1(isnan(x1))=0;
x2(isnan(x2))=0;

[ps1(cnt),dAUC1(cnt)] = logrank(x1,x2,'High drivers','Low drivers');
text(1000,0.2,['P = ',num2str(ps1(cnt))])


canc=DFIN1.canc;
OS = DFIN1.OS;
death = DFIN1.death;

[a,b,c] = intersect(DFIN1.gene,PANDRIVE);
DRC = sum(DFIN1.MUT(b,:));


ss1=[];

ss1 = DRC>THR;

x1 = [OS(ss1==1),death(ss1==1)];
x2 = [OS(ss1==0),death(ss1==0)];
x1(isnan(x1))=0;
x2(isnan(x2))=0;
subplot(2,4,cnt+4)

[ps2(cnt),dAUC2(cnt)] = logrank(x1,x2,'High drivers','Low drivers');
text(1000,0.2,['P = ',num2str(ps2(cnt))])
cnt=cnt+1
end