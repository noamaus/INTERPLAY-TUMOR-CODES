%%% Plot Figure 1D
load('DFIN.mat') %%PAN CANCER data
load('PANDRIVE.mat')

[a,b,c] = intersect(DFIN.gene,PANDRIVE);



uc=unique(DFIN.canc);

RR1=[];PP1=[];
for i = 1:length(uc)
    i
    cc=strcmp(DFIN.canc,uc(i));
    for j = 1:length(a)
        
        [RR1(j,i),PP1(j,i)]=corr(DFIN.AN(cc),DFIN.MUT(b(j),cc)','type','Spearman');
        
    end
end
  

T2=RR1';
T3=T2;
for j = 1:32
    T3(j,isnan(T3(j,:)))=mean(T3(j,~isnan(T3(j,:))));
end
cgObj2 = clustergram(T3','Cluster','row') %%%!!
R1 = str2num(char(cgObj2.ColumnLabels));
[sv,si] = sort(mean(T3),'descend');

load('mapGR.mat')
T4 = T3(R1,si);
h1=heatmap(a(si),uc(R1),T4)
h1.Colormap=mapGR
caxis([-0.2,0.2])