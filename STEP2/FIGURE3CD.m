%%%% Plot Figure 3C-D

load('DFIN1.mat') %%gastrointestinal and endometrial data
load('DFIN2.mat') %%all other tumors data
load('SIG1.mat') %%Repair signature
load('SIG2.mat') %%Apoptosis signature
load('DFIN.mat') %%PAN CANCER data

scc = (sum(DFIN.MUT(SIG2,:))+1)./(sum(DFIN.MUT(SIG1,:))+1);%%APOP/repair

[sv1,si1] = sort(sum(DFIN.MUT(SIG1,:)'),'descend');
[sv2,si2] = sort(sum(DFIN.MUT(SIG2,:)'));
SIG11=SIG1(si1);
SIG22=SIG2(si2);

%%Get top and bottom scores
SCC=scc([find(scc>2),find(scc<0.5)]);
AN2=DFIN.AN([find(scc>2),find(scc<0.5)]);
X11=DFIN.MUT([SIG11],[find(scc>2),find(scc<0.5)]);
X22=DFIN.MUT([SIG22],[find(scc>2),find(scc<0.5)]);


%%PLOT MAPS
subplot(1,9,1:3)
h1=heatmap(DFIN.gene(SIG11),1:length(X11),X11')
load('mapR.mat')
h1.Colormap=mapR

subplot(1,9,5:6)
h1=heatmap(DFIN.gene(SIG22),1:length(X22),X22')
load('mapB.mat')
h1.Colormap=mapB

subplot(1,9,8)
h1=heatmap(AN2)
h1.Colormap=mapR
caxis([0,39])

subplot(1,9,9)
h1=heatmap(SCC')
load('mapG.mat')
h1.Colormap=mapG
