clear all
clc
%% load original data
load('Y_nopreprocess.mat','Y')
Xorig=tensor(Y.data);

%% add noise 
eta=0.3; %the level of noise
%   N=tensor(randn(size(Xorig)));
load ('N.mat','N')
%   save ('N.mat','N')
Xnoise=Xorig+eta*N/norm(N)*norm(Xorig);
%% set missing values
mis_perc=0.2*ones(1,20);
%  Wall=cell(1,length(mis_perc));
load('Wall.mat','Wall')
for q=1:length(mis_perc) % loop for 20 samples 
    % W=create_missing_data_pattern(s,mis_perc(q));
    % Wall{1,q}=W;
    W=Wall{1,q};
    data.Wall{q}=W;
    % set missing value to be nan and then preprocessing
    Xmiss=Xnoise;
    Xmiss(find(W.data==0))=nan;
    %% preprocessing
    %  centering across the condition mode
    XX=Xmiss.data;
    temp = XX(:,:);
    temp_centered = temp - repmat(nanmean(temp),size(temp,1),1);
    XX_centered = reshape(temp_centered, size(XX));
    % scaling in the second mode - using std or root mean square
    X_centered=tensor(XX_centered);
    X=X_centered;
    for j=1:size(X,2)
        temp = squeeze(X.data(:,j,:));
        rms = sqrt(nanmean((temp(:).^2)));
        XX(:,j,:) = temp/rms;           
    end
    Xpre=tensor(XX);
    Xpre(find(W.data==0))=0;
    %% Perform the CP model
    X=Xpre;
    nb_starts =40;
    nm_comp=2;
    optionsCP.factr=1e-10;
    optionsCP.maxIts = 10000;
    optionsCP.maxTotalITs=50000;
    optionsCP.printEvery  = 10000;
    Low{1}=-Inf*ones(size(X,1),nm_comp);
    Low{2}=-Inf*ones(size(X,2),nm_comp);
    Low{3}=zeros(size(X,3),nm_comp);
    goodness_X1 = strings(nb_starts,1);
    goodness_X = zeros(nb_starts,2); %Stores ExitMsg, Fit, F(error for lbfgsb)
    Fac_X = cell(nb_starts,1);
    out_X = cell(nb_starts,1);
    
    %Call cp_opt to fit the CP model
    
    for i=1:nb_starts
        if i==1
           
            [Fac_X{i}, ~, out_X{i}] =cp_wopt(X,W,nm_comp,'init','nvecs','lower',Low,'opt_option',optionsCP);
        else
           
            [Fac_X{i}, ~, out_X{i}] =cp_wopt(X,W,nm_comp,'init','randn','lower',Low,'opt_option',optionsCP);
            
        end
        
        goodness_X1(i) = out_X{i}.ExitMsg;
        goodness_X(i,2) = out_X{i}.OptOut.err(end,1);
    end
    
    %%
    [ff, index] = sort(goodness_X(:,2),'ascend');
    Fac_sorted = Fac_X(index);
    out_sorted = out_X(index);
    Fac_X_best = Fac_sorted{1};
    out_X_best = out_sorted{1};
    Fac = Fac_X_best;
    data.fit(q)=100- (norm(tensor(full(Fac)-X))^2/norm(X)^2*100);
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
% CP2_Miss20_data_noise_03_noit=data;
% save ('CP2_Miss20_data_noise_03_noit','CP2_Miss20_data_noise_03_noit')

% CP2_Miss20_data_noise_01_noit=data;
% save ('CP2_Miss20_data_noise_01_noit','CP2_Miss20_data_noise_01_noit')
%%% save ('Wall.mat','Wall')

