% This is an example script for missing data estimation through cross-validation using the CP model.

% We use Tensor Toolbox as well as the L-BFGS-B implementation from https://github.com/stephenbeckr/L-BFGS-B-C
% In addition, parts of the scripts may require the dataset object (https://eigenvector.com/software/dataset-object/), publically available.

%% load data
load LOS_beta001_alpha05.mat
Xorig = tensor(Y.data);

%% add noise
eta    = 0.3; %the level of noise
N      = tensor(randn(size(Xorig)));
Xnoise = Xorig + eta*N/norm(N)*norm(Xorig);

%% algorithmic options
nb_starts = 20;
nm_comp   = 2;
optionsCP.maxIts      = 10000;
optionsCP.maxTotalITs = 50000;
optionsCP.printEvery  = 10000;
Low{1} = -Inf*ones(size(Xorig,1),nm_comp);
Low{2} = -Inf*ones(size(Xorig,2),nm_comp);
Low{3} = zeros(size(Xorig,3),nm_comp);

%% CP model with missing data
mis_perc = 0.2*ones(1,20); % missing percentage

for q = 1:length(mis_perc) % loop for 20 different missing data patterns
    
    W = create_missing_data_pattern(size(Xorig),mis_perc(q));
    data.Wall{q} = W; %store the missing position for each sample
        
    % set missing value to be NaN
    Xmiss = Xnoise;
    Xmiss(find(W.data==0))= NaN;
    
    %  preprocess the data
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
    Xpre = tensor(XX);
    Xpre(find(W.data==0))=0;
    
    %  Fit the CP model
    X = Xpre;
    goodness_X1 = strings(nb_starts,1);
    goodness_X  = zeros(nb_starts,1); %Stores ExitMsg, Fit, F(error for lbfgsb)
    Fac_X = cell(nb_starts,1);
    out_X = cell(nb_starts,1);
    
    % Call cp_wopt to fit the model
    for i = 1:nb_starts
        if i==1
            [Fac_X{i}, ~, out_X{i}] = cp_wopt(X,W,nm_comp,'init','nvecs','lower',Low,'opt_option',optionsCP);
        else
            [Fac_X{i}, ~, out_X{i}] = cp_wopt(X,W,nm_comp,'init','randn','lower',Low,'opt_option',optionsCP);
        end
        goodness_X1(i) = out_X{i}.ExitMsg;
        goodness_X(i)  = out_X{i}.OptOut.err(end,1);
    end
    
    % the factor for each sample
    [ff, index] = sort(goodness_X(:),'ascend');
    Fac = Fac_X{index(1)};
    data.exitMS{q} = goodness_X1(index(1));
    
    % To compute the TCS, we need the postprocessing
    % pull back the data
    Xtilde = full(Fac);
    X      = X_centered;
    % scaling back within the metabolites mode
    for j=1:size(X,2)
        temp      = squeeze(X(:,j,:));
        temptilde = squeeze(Xtilde.data(:,j,:));
        rms       = sqrt(nanmean((temp(:).^2)));
        XX(:,j,:) = temptilde*rms;
    end
    Xtilde = tensor(XX);
    % centering back across the subjects mode
    temptilde = Xtilde.data(:,:);
    X = Xmiss;
    temp = X.data(:,:);
    temptilde_uncentereed = temptilde+repmat(nanmean(temp),size(temp,1),1);
    Xtilde_uncenter       = reshape(temptilde_uncentereed, size(X));
    Xpredic = tensor(Xtilde_uncenter);
    
    data.TCS(q)  = norm(tensor(times(1-W,Xpredic-Xnoise)))/norm(tensor(times(1-W,Xnoise)));
    data.TS(q)   = norm(tensor(times(W,Xpredic-Xnoise)))/norm(tensor(times(W,Xnoise)));
    data.Xhat{q} = Xpredic;
end


%% plot actual vs. predicted for each missing data pattern
data.X = Xorig;
figure
for q = 1:length(data.Xhat)
    plot(data.X(find(data.Wall{q}==0)),data.Xhat{q}(find(data.Wall{q}==0)),'*','Markersize',11)
    set(gca,'fontsize',20)
    xlabel('actual')
    ylabel('predicted')
    hold on;
end

% %% save the data
% CP2_Miss20_data_noise_03=data;
% save ('CP2_Miss20_data_noise_03','CP2_Miss20_data_noise_03')

