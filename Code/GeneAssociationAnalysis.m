clc
clear
% load('GroupD_Value.mat')
load('MV_Harmonized_Zscore_D_value.txt')
load('parcelExpressionExcluUnAssign.mat')
%% PLS analysis

X=parcelExpressionExcluUnAssign(:,2:end);
Tem = parcelExpressionExcluUnAssign(:,1);
Y =MV_Harmonized_Zscore_D_value(Tem);
X=zscore(X);
Y=zscore(Y);
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y);

% plot explained ratio
dim=15;
plot(1:dim,100*PCTVAR(2,1:dim),'-o','LineWidth',2,'Color',[140/255,0,0]);
set(gca,'linewidth',2.7);
set( gca, 'Position', [ 0.005, 0.007, 0.98, 0.98 ] );
box off

% correlation between the first four pls component scores and Y
[R1,p1]=corr(XS(:,1),Y);
[R2,p2]=corr(XS(:,2),Y);
[R3,p3]=corr(XS(:,3),Y);
[R4,p4]=corr(XS(:,4),Y);
save('PLSscoreY_R.mat', 'R1','R2','R3','R4')
%% spatial autocorrelation corrected permutation test to assess the significance of PLS component variance explained ratios

load('surrogate_maps.mat')
for dim=1:15
   dim 
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
tem=cumsum(100*PCTVAR(2,1:dim));
Rsquared = tem(dim);

    for counter=1:size(surrogate_maps,1)
        counter
        Yp=surrogate_maps(counter,:)';
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,dim);
        tem=cumsum(100*PCTVAR(2,1:dim));
        Rsq(counter) = tem(dim);
    end

R(dim)=Rsquared;
p(dim)=length(find(Rsq>=Rsquared))/size(surrogate_maps,1);
end
save('PLS_VarianceExplained.mat','R','p')

%% get the ordered gene list

load('100DS512scaledRobustSigmoidNSGRNAseqQC1Lcortex_ROI_NOdistCorrEuclidean.mat')
genesSymbol = probeInformation.GeneSymbol;
geneID = probeInformation.EntrezID;

% Do PLS in 1 dimensions (with 1 components):
dim = 1;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
PLS1_Group_D = [XS(:,1),Y];
save('PLS1_Group_D.mat','PLS1_Group_D')

%align PLS components with desired direction for interpretability 
if R1<0  
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end

%Order the genes accroding to the initial PLS weights
[PLS1w,Indx1] = sort(stats.W(:,1),'descend');
PLS1w_GeneSymbol=genesSymbol(Indx1);
PLS1w_GeneID=geneID(Indx1);
PLS1weights=[];

bootnum=1000;
for counter=1:bootnum
    counter
    myresample = randsample(size(X,1),size(X,1),1); 
    Xr=X(myresample,:); 
    Yr=Y(myresample,:); 
  
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim);
    
    tem=stats.W(:,1);
%     order the newly obtained weights the same way as initial PLS 
    newWei=tem(Indx1); 
%     As the sign of PLS components is arbitrary, make sure this aligns between runs
    if corr(PLS1w,newWei)<0 
        newWei=-1*newWei;
    end
    PLS1weights=[PLS1weights,newWei];
    
end

%standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');

%z-score weights
temp1=PLS1w./PLS1sw';

%rank gene according to bootstrap weights
[Z1,ind1]=sort(temp1,'descend');
OrderedGeneSymbol=PLS1w_GeneSymbol(ind1);
OrderedGeneID=PLS1w_GeneID(ind1);
save('OrderedZvalue.mat','Z1');
save('OrderedGeneID.mat','OrderedGeneID');
save('OrderedGeneSymbol.mat','OrderedGeneSymbol');






