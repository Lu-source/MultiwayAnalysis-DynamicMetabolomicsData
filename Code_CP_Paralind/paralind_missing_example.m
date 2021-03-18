
%%
clear all
addpath('/Users/ll/Documents/MATLAB/toolbox/PLS_Toolbox_891')
save path 
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
Fac_X_best = Fac_X{index(1)};
Fac = Fac_X_best;
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
% Paralind2_Miss20_data_noise_03_noit=data;
% save ('Paralind2_Miss20_data_noise_03_noit','Paralind2_Miss20_data_noise_03_noit')