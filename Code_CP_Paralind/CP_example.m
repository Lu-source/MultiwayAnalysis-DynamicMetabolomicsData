% This script is an example script showing how to fit a CP model to the data generated 
% through simulations of dynamic systems. 

% We use Tensor Toolbox as well as the L-BFGS-B implementation from https://github.com/stephenbeckr/L-BFGS-B-C
% In addition, parts of the scripts may require the dataset object (https://eigenvector.com/software/dataset-object/), publically available.
% For core consistency computation, we also use the corcord function from the Nway toolbox (http://www.models.life.ku.dk/nwaytoolbox).

%% load data: 
% Y is a subjects by metabolites by time points tensor.
load('LOS_beta03_alpha05.mat','Y')
X = tensor(Y.data);
labelss = Y.label{2};
s=size(X);


%% preprocess data
%  centering across the subject mode
XX   = X.data;
temp = XX(:,:);
temp_centered = temp - repmat(mean(temp),size(temp,1),1);
XXX = reshape(temp_centered, size(XX));
% scaling in the metabolites mode - using root mean square
for j=1:size(X,2)
    temp = squeeze(XXX(:,j,:));
    rms  = sqrt(mean((temp(:).^2)));
    XX(:,j,:) = temp/rms;          
end
X    = tensor(XX);


%% plot the preprocessed data
figure
for i=1:s(2)
    subplot(4,3,i)
    for j=1:s(1)
        if j<=10
            plot(1:s(3),squeeze(X.data(j,i,:)),'r','Linewidth',1)
            xlim([1 20]);
            hold on
        elseif j<=20
            plot(1:s(3),squeeze(X.data(j,i,:)),'b','Linewidth',1)
            xlim([1 20]);
            hold on
        end
    end
    title(labelss(i,:))
    set(gca,'FontSize', 14)
end


%% CP model
nb_starts = 60;
nm_comp   = 2;
optionsCP.maxIts = 10000;
optionsCP.maxTotalITs = 50000;
optionsCP.printEvery  = 10000;
Low{1}=-Inf*ones(size(X,1),nm_comp);
Low{2}=-Inf*ones(size(X,2),nm_comp);
Low{3}=zeros(size(X,3),nm_comp);    % nonnegativity constraints in the third mode
goodness_X1 = strings(nb_starts,1); % Stores ExitMsg for CP model
goodness_X  = zeros(nb_starts,2);   % Stores  Fit and  F(error for lbfgsb)
Fac_X = cell(nb_starts,1);
out_X = cell(nb_starts,1);
for i=1:nb_starts
    if i==1
        [Fac_X{i}, ~, out_X{i}] = cp_opt(X,nm_comp,'init','nvecs','lower',Low,'opt_option',optionsCP);
    else
        [Fac_X{i}, ~, out_X{i}] = cp_opt(X,nm_comp,'init','randn','lower',Low,'opt_option',optionsCP);
    end
    goodness_X1(i)  = out_X{i}.ExitMsg;
    goodness_X(i,1) = out_X{i}.Fit;
    goodness_X(i,2) = out_X{i}.OptOut.err(end,1);
end

%% uniqueness test
% 0 -> NOT unique
% 1 -> Unique
% 2 -> Inconclusive, need more random starts
good_flag = find(goodness_X1(:) == 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL.' | goodness_X1(:) == 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH.');
if length(good_flag)>=1
    F_round      = round(goodness_X(good_flag,2),8);
    best_F_index = good_flag(F_round == min(F_round));
    if length(best_F_index) < 2
        F_round      = round(goodness_X(good_flag,2),5);
        best_F_index = good_flag(F_round == min(F_round));
    end   
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

%% fit, CC, TC, C12 values
uniqueness        = unique_test;
Fac = Fac_X{best_F_index(1)};
out_X_best = out_X{best_F_index(1)};
fit  = out_X_best.Fit; %model fit
Consistency = corcond(X.data,normalize(Fac,1),[],0); %core consistency
tc   = TC(Fac.U);      %tucker congruence
C12  =(Fac.U{1}(:,1)'*Fac.U{1}(:,2))/norm(Fac.U{1}(:,1))/norm(Fac.U{1}(:,2)); 


%%  plot CP model
Z2={'subjects','metabolites','time'};
Z3{1}=1:s(1);Z3{2}=1:s(2);Z3{3}=1:s(3);
labels{1} = Y.label{1};
labels{2} = Y.label{2};
labels{3} = Y.label{3};
Leglab = {'Comp1', 'Comp2'};

figure % scatter plot for the subjects mode
for i=1:3
    if i==1
        plot(Fac.U{i}(1:10,1),Fac.U{i}(1:10,2),'ro','Markersize',11,'MarkerFaceColor','r','LineWidth',1.5)
        hold on
        plot(Fac.U{i}(11:20,1),Fac.U{i}(11:20,2),'bo','Markersize',11,'MarkerFaceColor','b','LineWidth',1.5)
        text(Fac.U{i}(:,1),Fac.U{i}(:,2),labels{i},'VerticalAlignment','bottom','HorizontalAlignment','right','fontsize',15)
        set(gca,'FontSize', 20)
        xlabel('component 1','FontSize', 10)
        ylabel('component 2','FontSize', 10)
        title('subjects')
        
    end
    set(gca,'FontSize', 18)
    
end
figure % factor plot for the metabolites mode and time mode
for i=2:3
    subplot(2,1,i-1)
    if i==2
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




