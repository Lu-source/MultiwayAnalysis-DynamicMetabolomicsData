% This is an example script for missing data estimation through cross-validation using the Paralind model.

%To run the code, 'paralind_Lu_ortho.m' is needed.
% 'paralind_Lu_ortho.m' is an improved version of  the 'paralind.m' from
% http://www.models.life.ku.dk/paralind, with additional orthogonal
% constraint added to the factors and with bugs fixed for missing data


% In addition, parts of the scripts may require the dataset object
%(https://eigenvector.com/software/dataset-object/), publically available.


%% load data
load('GLM_beta002_PFKalpha05.mat','Y')
Xorig=tensor(Y.data);

%% add noise
eta=0.3; %the level of noise
N=tensor(randn(size(Xorig)));
Xnoise=Xorig+eta*N/norm(N)*norm(Xorig);

%% algorithmic options
nb_starts = 40;
H0 =[1 1];
A0 = [];B0 = [];C0 = [];
Options(1) = 1e-10;
Options(2) = 10000;
Options(3) = 2;
R=1;S=2;
Constraints = [0 -1 3 1];

%% Paralind model with missing data
mis_perc=0.2*ones(1,20); % missing percentage

for q=1:length(mis_perc) % loop for 20 different missing data patterns
    
    W=create_missing_data_pattern(size(Xorig),mis_perc(q));
    data.Wall{q}=W; %store the missing position for each sample
    
    % set missing value to be nan
    Xmiss=Xnoise;
    Xmiss(find(W.data==0))=nan;
    
    % preprocess the data
    %  centering across the subjects mode
    XX   = Xmiss.data;
    temp = XX(:,:);
    temp_centered = temp - repmat(nanmean(temp),size(temp,1),1);
    X_centered   = reshape(temp_centered, size(XX));
    %  scaling in the metabolites mode - using root mean square
    X = X_centered;
    for j = 1:size(X,2)
        temp = squeeze(X(:,j,:));
        rms  = sqrt(nanmean((temp(:).^2)));
        XX(:,j,:) = temp/rms;
    end
    X=tensor(XX);
    
    
    % perform the Paralind model
    FactorsXL=cell(1,nb_starts);
    Fac_X=cell(1,nb_starts);
    for i=1:nb_starts
        [A,H,B,C,fit,it,explainvar]=paralind_Lu_ortho(X.data,R,S,Constraints,Options,H0,A0,B0,C0);
        FactorsXL{i}{1}=A*H;
        FactorsXL{i}{2}=B;
        FactorsXL{i}{3}=C;
        Fac_X{i}=ktensor(FactorsXL{i});
        erF(i)=norm(times(W,X-full(Fac_X{i})));
        
        
    end
    
    % % the factor for each sample
    [ff, index] = sort(erF,'ascend');
    Fac = Fac_X{index(1)};
    
    % To compute the TCS, we need the postprocessing
    % pull back the data
    Xtilde=full(Fac);
    X=X_centered;
    %%scaling back within the metabolites mode
    for j=1:size(X,2)
        temp = squeeze(X(:,j,:));
        temptilde = squeeze(Xtilde.data(:,j,:));
        rms = sqrt(nanmean((temp(:).^2)));
        XX(:,j,:) = temptilde*rms;
    end
    Xtilde=tensor(XX);
    %%centering back across the condition mode
    temptilde=Xtilde.data(:,:);
    X=Xmiss;
    temp = X.data(:,:);
    temptilde_uncentereed=temptilde+repmat(nanmean(temp),size(temp,1),1);
    Xtilde_uncenter=reshape(temptilde_uncentereed, size(X));
    Xpredic=tensor(Xtilde_uncenter);
    
    data.TCS(q)=norm(tensor(times(1-W,Xpredic-Xnoise)))/norm(tensor(times(1-W,Xnoise)));
    data.TS(q)=norm(tensor(times(W,Xpredic-Xnoise)))/norm(tensor(times(W,Xnoise)));
    data.Xhat{q}=Xpredic;
end

%% plot predicted--actual
data.X=Xorig;
figure
for q=1:length(data.Xhat)
    plot(data.X(find(data.Wall{q}==0)),data.Xhat{q}(find(data.Wall{q}==0)),'*','Markersize',11)
    hold on
    set(gca,'fontsize',20)
    xlabel('actual')
    ylabel('predicted')
end
% %% save the data
% Paralind2_Miss20_data_noise_03=data;
% save ('Paralind2_Miss20_data_noise_03','Paralind2_Miss20_data_noise_03')
