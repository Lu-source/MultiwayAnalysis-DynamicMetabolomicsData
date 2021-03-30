
% An example script for computing the tensor completion score (TCS)
% using Paralind model when we perform cross-validation.
% To run this code, PLS toolbox is needed

%% load data
load('GLM_beta002_PFKalpha05.mat','Y')
Xorig=tensor(Y.data);

%% add noise
eta=0.3; %the level of noise
N=tensor(randn(size(Xorig)));
Xnoise=Xorig+eta*N/norm(N)*norm(Xorig);

%% Paralind model with missing data
mis_perc=0.2*ones(1,20); % missing percentage

for q=1:length(mis_perc) % loop for 20 samples
    
    W=create_missing_data_pattern(s,mis_perc(q));
    data.Wall{q}=W; %store the missing position for each sample
    
    %% set missing value to be nan
    Xmiss=Xnoise;
    Xmiss(find(W.data==0))=nan;
    
    %% preprocess the data 
    %  centering across the subjects mode
    XX=Xmiss.data;
    temp = XX(:,:);
    temp_centered = temp - repmat(nanmean(temp),size(temp,1),1);
    XX_centered = reshape(temp_centered, size(XX));
    X_centered=tensor(XX_centered);
    % scaling in the metabolites mode - using root mean square
    X=X_centered;
    for j=1:size(X,2)
        temp = squeeze(X.data(:,j,:));
        rms = sqrt(nanmean((temp(:).^2)));
        XX(:,j,:) = temp/rms;
    end
    Xpre=tensor(XX);
   
    
%% perform the Paralind model
X=Xpre;
nb_starts=4;
opt = parafac('options');
opt.init =10;
opt.plots = 'off';
opt.stopcrit(1)=1e-10;
opt.constraints{1}.type='rightprod';
opt.constraints{1}.advanced.linearconstraints.matrix = [1 1];
opt.constraints{3}.type='nonnegativity';
opt.constraints{2}.type='orthogonality';
FactorsXL=cell(1,nb_starts);
Fac_X=cell(1,nb_starts);
nm_comp=2;
for i=1:nb_starts
    
    m{i} = parafac(X.data, nm_comp, opt);
    FactorsXL{i}{1}=m{i}.loads{1};
    FactorsXL{i}{2}=m{i}.loads{2};
    FactorsXL{i}{3}=m{i}.loads{3};
    Fac_X{i}=ktensor(FactorsXL{i});
     erF(i)=norm(X-full(Fac_X{i})); 
end
%%
[ff, index] = sort(erF,'ascend');
Fac = Fac_X{index(1)};

%% To compute the TCS, we need the postprocessing
    % pull back the data
    Xtilde=full(Fac);
    X=X_centered;
    %%scaling back within the metabolites mode
    for j=1:size(X,2)
        temp = squeeze(X.data(:,j,:));
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
    Xhat=Xpredic;
    data.Xhat{q}=Xhat;
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
%% save the data
% Paralind2_Miss20_data_noise_03=data;
% save ('Paralind2_Miss20_data_noise_03','Paralind2_Miss20_data_noise_03')