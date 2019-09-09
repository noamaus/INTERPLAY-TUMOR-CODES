%%% Plot Figure 1C

load('PANDRIVE.mat')
load('DFIN1.mat') %%gastrointestinal and endometrian data
load('DFIN2.mat') %%all other tumors data

cnt=1;
ps1=[];ps2=[];dAUC1=[];dAUC2=[];

[a,b,c] = intersect(DFIN1.gene,PANDRIVE);



%%% For each threshold 1-15, get KM P-value and dAUC based on #drivers, for
%%% ech tumor cluster
for i = 1:15

THR=i;

canc=DFIN2.canc;
OS = DFIN2.OS;
death = DFIN2.death;
DRC = sum(DFIN2.MUT(b,:));
ss1 = DRC>THR;


x1 = [OS(ss1==1),death(ss1==1)];
x2 = [OS(ss1==0),death(ss1==0)];
x1(isnan(x1))=0;
x2(isnan(x2))=0;

[ps1(cnt),dAUC1(cnt)] = logrank2(x1,x2,'High drivers','Low drivers');


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

[ps2(cnt),dAUC2(cnt)] = logrank2(x1,x2,'High drivers','Low drivers');%%Gastrointestinal,endometrial
cnt=cnt+1
end


%% PLOT the dAUC for each cluster
ps1(ps1==0) = 1e-20;

stem([dAUC1])
hold on 
stem([dAUC2])
hold on
ylabel('Kaplan-Meier dAUC')
scatter(1:15,dAUC2,-log(ps2)*10,'filled')
scatter(1:15,dAUC1,-log(ps1)*10,'filled')

legend('Other','Gastrointestinal,endometrial')



