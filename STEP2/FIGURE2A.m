%%% Plot Figure 2A
load('DDRGENES.mat') %% List of all DDR genes
load('DFIN.mat') %%PAN CANCER data
load('RESPL.mat') %%% DD Selection scores
load('SIG1.mat') %%Repair signature
load('SIG2.mat') %%Apoptosis signature


uc=unique(DFIN.canc)
[aa,bb,cc] = intersect(DFIN.gene(SIG1),DFIN.gene(DDRGENES));
[aa2,bb2,cc2] = intersect(DFIN.gene(SIG2),DFIN.gene(DDRGENES));


id1=[30,6,26,23,31]; %% gastrointestinal, endometrial indices
id2=[1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,27,28,29,32]; %% other tumors indices


T11 = -log([RESPL.FINPV1(cc,id1)]);
T12 = -log([RESPL.FINPV1(cc2,id1)]);

T21 = -log([RESPL.FINPV2(cc,id2)]);
T22 = -log([RESPL.FINPV2(cc2,id2)]);

[svv1,sii1] = sort(mean(T11'),'descend');
[svv2,sii2] = sort(mean(T22'),'descend');

[sc1,sic1] = sort(mean(T11),'descend');
[sc2,sic2] = sort(mean(T22),'descend');

T11 = -log([RESPL.FINPV1(cc(sii1),id1(sic1))]);
T12 = -log([RESPL.FINPV1(cc2(sii2),id1(sic1))]);

T21 = -log([RESPL.FINPV2(cc(sii1),id2(sic2))]);
T22 = -log([RESPL.FINPV2(cc2(sii2),id2(sic2))]);


load('map2.mat')
subplot(2,2,1)
h1=heatmap(uc(id1(sic1)),[aa(sii1)],T11)
h1.Colormap=map2
caxis([-0.5,4.5])

subplot(2,2,2)
h1=heatmap(uc(id1(sic1)),[aa2(sii2)],T12)
h1.Colormap=map2
caxis([-0.5,4.5])

subplot(2,2,3)
h1=heatmap(uc(id2(sic2)),[aa(sii1)],T21)
h1.Colormap=map2
caxis([-0.5,4.5])

subplot(2,2,4)
h1=heatmap(uc(id2(sic2)),[aa2(sii2)],T22)
h1.Colormap=map2
caxis([-0.5,4.5])