
% An example script showing how to fit a Paralind model to the data generated
% by dynamic metabolomics data sets.

%% load data
load('GLM_beta002_PFKalpha05.mat','Y')
X=tensor(Y.data);
labelss = Y.label{2};
s=size(X);

%% preprocess data
%  centering across the subject mode
XX=X.data;
temp = XX(:,:);
temp_centered = temp - repmat(mean(temp),size(temp,1),1);
XXX = reshape(temp_centered, size(XX));
X=tensor(XXX);
% scaling in the metabolites mode - using root mean square
for j=1:size(X,2)
    temp = squeeze(X.data(:,j,:));
    rms = sqrt(mean((temp(:).^2)));
    XX(:,j,:) = temp/rms;
end
Xpre=tensor(XX);
X=tensor(Xpre);

%% plot the preprocessed data
figure
for i=1:s(2)
    subplot(4,3,i)
    for j=1:s(1)
        if j<=10
            plot(1:1:20,squeeze(X.data(j,i,:)),'r','Linewidth',1)
            xlim([1 20]);
            hold on
        elseif j<=20
            plot(1:1:20,squeeze(X.data(j,i,:)),'b','Linewidth',1)
            xlim([1 20]);
            hold on
        end
    end
    title(labelss(i,:))
    set(gca,'FontSize', 14)
end

%% Paralind model
nb_starts =60;
H0=[1 1];
A0=[];B0=[];C0=[];
Options(1)=1e-10;
Options(2)=10000;
Options(3)=2;
R=1;S=2;
Constraints = [0 -1 3 1];
for i=1:nb_starts
    [A,H,B,C, ~, ~,explainvar]=paralind_Lu_ortho(X.data,R,S,Constraints,Options,H0,A0,B0,C0);
    FactorsXL{i}{1}=A*H;
    FactorsXL{i}{2}=B;
    FactorsXL{i}{3}=C;
    Fac_X{i}=ktensor(FactorsXL{i});
    erF(i)=norm(X-full(Fac_X{i}));
    expvar(i)=explainvar;
end
%% uniqueness test
Fit_round = round(erF,8);
best_F_index = find(Fit_round == min(Fit_round));
if length(best_F_index) < 2
    F_round = round(goodness_X(good_flag,2),5); % can be adjusted if 5 is too tight
    best_F_index = good_flag(F_round == min(F_round));
end
eps = .05;
if length(best_F_index)==1
    unique_test = 2;
    disp('Need more random starts to determine uniqueness')
    worst_check = 0;
elseif length(best_F_index) > 1
    check_matrix = zeros(length(best_F_index));
    for i = 1:length(best_F_index)
        for j = 1:length(best_F_index)
            check_matrix(i,j) = score(Fac_X{best_F_index(j)},Fac_X{best_F_index(i)},'lambda_penalty',false);
        end
    end
    worst_check = min(min(check_matrix));
    if worst_check < (1-eps)
        unique_test = 0;
    else
        unique_test = 1;
    end
end


%%
unique=unique_test
[er,best_F_index]=sort(erF,'ascend');
Fac=Fac_X{best_F_index(1)};
fit=expvar(best_F_index(1))
% normalise the factor
nm_comp=S; % number of components
for j=1:nm_comp
    Fac.lambda(j)=1;
    for i=1:length(s)
        a=Fac.U{i}(:,j);
        Fac.lambda(j)=Fac.lambda(j)*Fac.lambda(j)*norm(a);
        Fac.U{i}(:,j)=a/norm(a);
    end
end
%%  plot CP model
Z2={'subjects','metabolites','time'};
Z3{1}=1:s(1);Z3{2}=1:s(2);Z3{3}=1:s(3);
labels{1} = Y.label{1};
labels{2} = Y.label{2};
labels{3} = Y.label{3};
Leglab = {'Comp1', 'Comp2'};

figure 
for i=1:3
    subplot(3,1,i)
    if i==1
        for j=1:nm_comp
            plot(Z3{i},Fac.U{i}(:,j),'-o','LineWidth',2.4)
            hold on
            title(Z2{i})
        end
    elseif i==2
        for j=1:nm_comp
            plot(Z3{i},Fac.U{i}(:,j),'-o','LineWidth',2.4)
            set(gca,'xtick',Z3{i},'xticklabel',labelss)
            xtickangle(90)
            hold on
            ax = gca;
            ax.YGrid = 'on';
            ax.GridAlpha = 1;
            title(Z2{i})
        end
    else
        for j=1:nm_comp
            plot(Z3{i},Fac.U{i}(:,j),'-o','LineWidth',2.4)
            hold on
            title(Z2{i})
        end
    end
    set(gca,'FontSize', 18)
    legend(Leglab,'TextColor','blue')
    
end

