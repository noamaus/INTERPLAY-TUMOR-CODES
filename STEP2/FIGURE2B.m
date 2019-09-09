load('DFIN1.mat') %%gastrointestinal and endometrial data
load('DFIN2.mat') %%all other tumors data
load('SIG1.mat') %%Repair signature
load('SIG2.mat') %%Apoptosis signature
load('DFIN.mat') %%PAN CANCER data




[sv1,si1] = sort(sum(DFIN.MUT(SIG1,:)'),'descend');
[sv2,si2] = sort(sum(DFIN.MUT(SIG2,:)'));


sc1 = (sum(DFIN1.MUT(SIG1,:)));%% Repair load, gastrointestinal and endometrial
sc2 = (sum(DFIN2.MUT(SIG2,:)));%% Apoptosis load, other tumors
load('PANDRIVE.mat')
[a,b,c] = intersect(DFIN.gene,PANDRIVE);
DRC1 = sum(DFIN1.MUT(b,:));%% Driver load, gastrointestinal and endometrial
DRC2 = sum(DFIN2.MUT(b,:));%% Driver load, other tumors


%%% Get correlation between signatures load and drivers loads in each
%%% cluster
[RC11,PC11] = corr(sc1',DRC1','type','Spearman');
[RC22,PC22] = corr(sc2',DRC2','type','Spearman');


subplot(2,1,1)
scatter(sc1',DRC1',70,'filled')
xlabel('Signature 1')
ylabel('Drivers load')
legend(['rho = ',num2str(RC11),'P = ',num2str(PC11)])

subplot(2,1,2)
scatter(sc2',DRC2',70,'filled')
legend(['rho = ',num2str(RC22),'P = ',num2str(PC22)])
xlabel('Signature 2')
ylabel('Drivers load')
