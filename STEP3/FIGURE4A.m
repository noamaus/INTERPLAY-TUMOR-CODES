%%%% Plots Figure 4A

load('DFIN11.mat') %% load gastrointestinal and endometrial data with MSI
load('SIG3.mat') %%  MSI-Aneuploidy signature



uc=unique(DFIN11.canc);
sm=sum(DFIN11.MUT(SIG3,:));  %% MSI-Aneuploidy signature load
sm2=sum(DFIN11.MUT(SIG3,:)'); %% MSI-Aneuploidy sum of mutation

[sv,si] = sort(sm,'descend');
[sv2,si2] = sort(sm2,'descend');


subplot(4,1,1:2)
X = DFIN11.MUT(SIG3(si2),si);
h1 = heatmap(DFIN11.sample(si),DFIN11.gene(SIG3(si2)),X)


subplot(4,1,3)
MSI2=strcmp(DFIN11.MSI,'MSI-H');
h2 = heatmap(double(MSI2(si))')
load('mapG.mat')
h2.Colormap=mapG


subplot(4,1,4)
h3 = heatmap(DFIN11.AN(si)')
load('mapR.mat')
h3.Colormap=mapR
caxis([1,39])



