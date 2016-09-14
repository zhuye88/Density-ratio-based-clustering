function [ F,match] = FmeasOPTICS( SetOfClusters,order,class )
%FMEASOPTICS calculate F measure for OPTICS results
%   Input  : X  : data set (m,n); m-objects, n-variables 
%            threshold :  density-ratio threshold for DBSCAN
%            Eps :  Eps neigbourhood for ReCon-DBSCAN
%            Eta :  Eta neigbourhood for ReCon-DBSCAN    
%            Matrix :  dissmilarity matrix
%   Output : class: cluster labels (m by 1)
%            type :  label types (1: core; 2: boundary; 3: noise)
 %% confusion matrix
matrix=zeros(max(class),length(SetOfClusters));


for i=1:length(SetOfClusters)
    Nclass=class-class;
    Tclass=Nclass;
    Nclass(SetOfClusters(i).start:SetOfClusters(i).end)=i;
    Tclass(order)=Nclass;
    index=find(Tclass==i);
    cindex=class(index);
    for j=1:length(index)
        matrix(cindex(j),i)=matrix(cindex(j),i)+1;
    end
end

%% recall
recall=matrix;
for j = 1:size(matrix,1)
        recall(j,:) = recall(j,:)/(sum(class==j)+0.0000001); % calculate the total positive examples 
end

%% precision
precision=matrix;
sumcol=sum(matrix);
if size(matrix,1)==1
    sumcol=matrix;
end

for j = 1:size(matrix,2)
        precision(:,j) = precision(:,j)/(sumcol(j)+0.0000001); % calculate the total positive examples 
end
%% fmeasure
fmeasure=2*precision.*recall./(precision+recall+0.0000001);

[match, cost] = hungarian(-fmeasure);

fmeasure=fmeasure(match==1);
F=sum(fmeasure)/max(class);
end

