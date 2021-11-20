function [Corr,FeatureIndex,PredictedScore]=SVR_LOOCV_Prediction( X_Features,Y_Scores,pthresh,Covariates,OutputPath)
% Leave one out cross validation 

cd 'MyToolbox/libsvm-3.24/matlab'
 FeatureIndex = cell(length(Y_Scores),1);
 PredictedScore = zeros(length(Y_Scores),1);
for counter = 1:length(Y_Scores)
 
 TestData = X_Features(counter,:);
 TestScore = Y_Scores(counter,1);
 TrainData = X_Features;
 TrainScores = Y_Scores;
 TrainData(counter,:)=[];
 TrainScores(counter)=[];
 %% regress covariates
 
 % prepare covariates
 Covariates_Test = Covariates(counter,:);
 Covariates_Train = Covariates;
 Covariates_Train(counter,:)=[];
 
 % demean covariates
 Covariates_Mean = mean(Covariates_Train,1);
 Covariates_Train = Covariates_Train-repmat(Covariates_Mean,length(Covariates_Train),1);
 Covariates_Test = Covariates_Test-Covariates_Mean;
 
 % regress covariates from behavior
 b = regress(TrainScores,[ones(length(TrainScores),1), Covariates_Train]);
 TrainScores = TrainScores-Covariates_Train*b(2:end);
 TestScore = TestScore-Covariates_Test*b(2:end);
 
 % regress covariates from features
 for counter_2 = 1:size(TrainData,2)
     b = regress(TrainData(:,counter_2),[ones(size(TrainData,1),1) Covariates_Train]);
     TrainData(:,counter_2) = TrainData(:,counter_2)-Covariates_Train*b(2:end);
     TestData(:,counter_2) = TestData(:,counter_2)-Covariates_Test*b(2:end);
 end
 
 % scale data
 MeanValue = mean(TrainData);
 STD = std(TrainData);
 TrainData = (TrainData-repmat(MeanValue,size(TrainData,1),1))./repmat(STD,size(TrainData,1),1);
 TestData = (TestData-MeanValue)./STD;
 
 % feature selection
 [~,p] = corr(TrainData,TrainScores);
 FeatureIndex{counter}= find(p<pthresh);
 TrainData = TrainData(:,FeatureIndex{counter});
 TestData = TestData( FeatureIndex{counter});
 
 % SVR training
 model = svmtrain(TrainScores,TrainData,'-s 3 -t 0 ');
 
 % SVR prediction
 [PredictedScore(counter,1),~,~] = svmpredict(TestScore,TestData,model);
 
end
Corr = corr(PredictedScore, Y_Scores);
save([OutputPath,filesep,'PredictedResultsRegressCov.mat'],'Corr','indx','FeatureWeight','PredictedScore')
end



