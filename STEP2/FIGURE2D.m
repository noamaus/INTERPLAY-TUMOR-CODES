%%% Plot FIGURE 2D

load('DFIN1.mat') %%gastrointestinal and endometrial data
load('DFIN2.mat') %%all other tumors data

load('SIG1.mat') %%Repair signature
load('SIG2.mat') %%Apoptosis signature


subplot(1,2,1)
%%% Get  cancer ratio of apotosis to repair gastrointestinal and
%%% endometrial tumors
scc = (sum(DFIN1.MUT(SIG2,:))+1)./(sum(DFIN1.MUT(SIG1,:))+1);%%APOP/repair

canc=DFIN1.canc;
OS = DFIN1.OS;
death = DFIN1.death;
AN = DFIN1.AN;

x1 = [OS(scc>1),death(scc>1)];%%%APOP>repair
x2 = [OS(scc<1),death(scc<1)];%%%APOP<repair

x1(isnan(x1))=0;
x2(isnan(x2))=0;

[ppp] = logrank(x1,x2,'Apoptosis','Repair');
text(1000,0.2,['P = ',num2str(ppp)])
text(8000,0.1,['n = ',num2str(length(x1))])
text(8000,0.5,['n = ',num2str(length(x2))])



subplot(1,2,2)
%%% Get  cancer ratio of apotosis to repair other tumors
scc = (sum(DFIN2.MUT(SIG2,:))+1)./(sum(DFIN2.MUT(SIG1,:))+1);%%APOP/repair

canc=DFIN2.canc;
OS = DFIN2.OS;
death = DFIN2.death;
AN = DFIN2.AN;

x1 = [OS(scc>1),death(scc>1)];%%%APOP>repair
x2 = [OS(scc<1),death(scc<1)];%%%APOP<repair

x1(isnan(x1))=0;
x2(isnan(x2))=0;

[ppp] = logrank(x1,x2,'Apoptosis','Repair');
text(1000,0.2,['P = ',num2str(ppp)])
text(8000,0.1,['n = ',num2str(length(x1))])
text(8000,0.5,['n = ',num2str(length(x2))])




