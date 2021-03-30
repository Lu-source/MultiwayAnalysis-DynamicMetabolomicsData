function [A,H,B,C,fit,it,explainvar] = paralind(X,R,S,Constraints,Options,H,A,B,C);

%     ___________________________________________________
%
%                  THE PARALIND MODEL
%     ___________________________________________________
% 
% Algorithm to fit the PARALIND model which is an advanced variant of the 
% normal PARAFAC model. It handles certain types of interactions between
% the components as typically occur in rank-deficient systems.
% See Bro, Harshman & Sidiropoulos, Journal of Chemometrics, 1999 for 
% details on application and Bro, Multi-way analysis in the Food Industry, 
% Ph.D. thesis, University of Amsterdam, 1998 for the first applications
% 
%
% The PARALIND model is given
% 
% Xk = A*H*Dk*B' + Ek, k = 1, .., K
% 
% Xk is a slab of data (I x J). K is the number of slabs. 
% A (I x R) are the scores or first-mode loadings. Dk is a diagonal 
% matrix that holds the k'th row of C (K x S) in its diagonal. C 
% (K x S) is the third mode loadings, B is a J x S loading matrix
% of the second mode, and H is an R x S interaction matrix
% defining the interactions between the R first mode loadings and
% the S loadings in B and C.
% 
% INPUT
% 
% X
%   Three-way that holds the data.
%      X(:,:,1) = X1; X(:,:,2) = X2; etc.  
%   If you have your data in an 'unfolded' two-way array of size
%   I x JK (the three-way array is I x J x K), then simply type
%   X = reshape(X,[I J K]); to convert it to an array.
%
% R
%   The number of components to extract in A
% S
%   The number of components to extract in B and C
% 
% Constraints
%   Vector of length 4. The first element defines constraints 
%   imposed in A, the second in H, third in B and fourth in C.
% 
%   If Constraints = [a b c d], the following holds. If 
%   a = -1 => keep A fixed by its initial values throughout
%   a =  0 => no constraints in A
%   a =  1 => nonnegativity in A
%   a =  2 => unimodality in A 
%   same holds for b, c, and d for H, B an C respectively (unimodal not for H)
%   a =  3 => orthogonality in A % by Lu, same holds for b for B
%
% Options
%   An optional vector of length 3
%   Options(1) Convergence criterion
%            1e-7 if not given or given as zero
%   Options(2) Maximal iterations
%            default 2000 if not given or given as zero
%   Options(3) Initialization method
%            A rather slow initialization method is used per default
%            but it pays to investigate in avoiding local minima.
%            Experience may point to faster methods (set Options(3)
%            to 1 or 2). You can also change the number of refits etc.
%            in the beginning of the m-file
%            0 => best of 20 runs of maximally 100 iterations (default)
%            1 => based on SVD
%            2 => random numbers
%            3 => input values
%   Options(4) Cross-validation
%            0 => no cross-validation
%            1 => cross-validation splitting in 7 segments
%   Options(5) show output
%            0 => show standard output on screen
%            1 => hide all output to screen
%
% AUXILIARY
% - Missing elements: Use NaN for missing elements
% - You can input initial values by using the input argument
%           (X,F,Constraints,Options,A,H,C,P);
%
% I/O
% 
% Short 
% [A,H,B,C]=paralind(X,R,S);
%
% Long
% [A,H,B,C,fit]=paralind(X,R,S,Constraints,Options);
%
% Copyright
% Rasmus Bro
% DK, 1998
% rb@life.ku.dk
%
% Reference to algorithm
% Bro, Harshman, Sidiropoulos, A new model for multi-way rank-deficient data,
% Journal of Chemometrics, 2007, Submitted

% $ Version 1.00 $ Date 19. April 1999 $ Not compiled $ 
% 
% This M-file and the code in it belongs to the holder of the 
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. Furthermore, the
% code can not be made part of any toolbox or similar.
% In case of doubt, contact the holder of the copyrights. 
%
% Rasmus Bro
% Chemometrics Group, Food Technology
% Department of Food and Dairy Science
% Royal Veterinary and Agricultutal University
% Rolighedsvej 30, DK-1958 Frederiksberg, Denmark
% Phone  +45 35283296
% Fax    +45 35283245
% E-mail rb@kvl.dk



dbg=0;

% NB Remove equality and fixed variables

if nargin==0    
    disp(' '),    
    disp(' '),
    disp(' THE THREE-WAY PARALIND MODEL'),    
    disp(' '),    
    disp(' Type <<help paralind>> for more info'),    
    disp('  '),
    disp(' I/O '),
   disp(' [A,H,B,C]=paralind(X,R,S);'),
   disp(' '),    
   disp(' Or optionally'),
   disp(' '),    
   disp(' [A,H,B,C,fit]=paralind(X,R,S,Constraints,Options);'),    
   disp(' '),    
   disp(' Options=[Crit MaxIt Init Xval Show]'),    
   disp(' '),    
   disp(' '),    
   return, 
elseif nargin<3&~all(X=='demo')    
    error(' The inputs X and R,S must be given') 
end 
clf

if isstr(X) & all(X=='demo')
   error('not implemented yet')
end


ShowFit    = 1000; % Show fit every 'ShowFit' iteration
delta      = 10; % For line-search
NumRep     = 20; %Number of repetead initial analyses
NumItInRep = 100; % Number of iterations in each initial fit
if ~length(size(X))==3
   error(' X must be a three-way array')
end

K = size(X,3);
%reset random number generators
rand('state',sum(100*clock))
randn('state',sum(100*clock))

if nargin<4
   Constraints = [0 0 0 0];
end
if length(Constraints)<4
   Constraints = [0 0 0 0];
   disp(' Length of Constraints must be four. It has been set to zeros')
end

if nargin < 5
   Options = zeros(1,5);
end
if length(Options)<5
   Options = Options(:);
   Options = [Options;zeros(5-length(Options),1)];
end

% Convergence criterion
if Options(1)==0
   ConvCrit = 1e-7;
else
   ConvCrit = Options(1);
end
if Options(5)==0
   disp(' ')
   disp(' ')
   disp([' Convergence criterion        : ',num2str(ConvCrit)])
end


% Maximal number of iterations 
if Options(2)==0
   MaxIt = 2000;
else
   MaxIt = Options(2);
end
ConstraintOptions=[ ...
      'Fixed                     ';...
      'Unconstrained             ';...
      'Non-negativity constrained';...
      'Unimodality constrained   ';...
      'Orthogonal constrained    '];
if Options(5)==0
   disp([' Maximal number of iterations : ',num2str(MaxIt)])
   disp([' Number of A factors          : ',num2str(R)])
   disp([' Number of B/C factors        : ',num2str(S)])
   disp([' Interaction matrix, H        : ',ConstraintOptions(Constraints(2)+2,:)])
   disp([' Loading 1. mode, A           : ',ConstraintOptions(Constraints(1)+2,:)])
   disp([' Loading 2. mode, B           : ',ConstraintOptions(Constraints(3)+2,:)])
   disp([' Loading 3. mode, C           : ',ConstraintOptions(Constraints(4)+2,:)])
   if length(Options)>9
      disp([' Equality in 3. mode, C       : ',num2str(Options(10))])
   end
   if length(Constraints)>5
      disp([' Zero-fixed parameters C      : '])
   end
   disp(' ')
end

% Initialization method
initi = Options(3);

% Fix first and last ten elements in C as described in Bro 98, thesis
if length(Constraints)>4
   Cindi = ones(K,S);
   Cindi(1:Constraints(5),[1 3 5])=0;
   Cindi(end-Constraints(5)+1:end,[2 4 6])=0;
else
   Cindi = ones(K,S);
end   

% Find missing and replace with average 
MissingElements = 0;
MissNum=0;AllNum=0;
MissingOnes=zeros(size(X));
for k = 1:K
   x=X(:,:,k);
   miss = sparse(isnan(x));
   MissingOnes(:,:,k) = miss;
   if any(miss(:))
      MissingElements = 1;
      % Replace missing with mean over slab (not optimal but what the heck)
      % Iteratively they'll be replaced with model estimates
      x(find(miss)) = mean(x(find(~miss)));
      X(:,:,k) = x;
      MissNum = MissNum + prod(size(find(miss)));
      AllNum = AllNum + prod(size(x));
   end
end
if MissingElements
   if Options(5)==0
      PercMiss = 100*MissNum/AllNum;
      RoundedOf = .1*round(PercMiss*10);
      disp([' Missing data handled by EM   : ',num2str(RoundedOf),'%'])
   end
end
clear x

% INTIALIZE LOADINGS
% Set non-existing loadings to []
if nargin<6,
    H=[];
end,
if nargin<7,
    A=[];
end,
if nargin<8,
    B=[];
end,
if nargin<9,
    C=[];
end
% Give initial ones unique name
Ho=H;Ao=A;Bo=B;Co=C;
if initi==0
   if Options(5)==0
      disp([' Using best of ',num2str(NumRep)]) 
      disp(' initially fitted models')
   end
   Opt = Options;
   Opt(1) = Options(1)/20;
   Opt(2) = NumItInRep; % Max NumItInRep iterations
   Opt(3) = 1;  % Init with SVD
   Opt(4) = 0;
   Opt(5) = 1;
   [A,H,B,C,bestfit] = paralind(X,R,S,Constraints,Opt,Ho,Ao,Bo,Co);
   AllFit = bestfit;
   for i = 2:NumRep
      Opt(3) = 2;   % Init with random
      [a,h,b,c,fit] = paralind(X,R,S,Constraints,Opt,Ho,Ao,Bo,Co);
      AllFit = [AllFit fit];
      if fit<bestfit
         A=a;H=h;C=c;B=b;
         bestfit = fit;
      end
   end
   AddiOutput.AllFit = AllFit;
   if Options(5)==0
      for ii=1:length(AllFit)
         disp([' Initial Model Fit            : ',num2str(AllFit(ii))])
      end
   end
   % Initialize by SVD
elseif initi==1
   if Options(5)==0
      disp(' SVD based initialization')
   end
   XtX = X(:,:,1)*X(:,:,1)';
   XXt = X(:,:,1)'*X(:,:,1);
   for k = 2:K
      XtX = XtX + X(:,:,k)*X(:,:,k)';
      XXt = XXt + X(:,:,k)'*X(:,:,k);
   end
   if all(size(A)~=[size(X,1) R]) % Generate A if not given beforehand
      [A,s,v]=svd(XtX,0);
      A=A(:,1:R);
   end
   if all(size(B)~=[size(X,2) S]) % Generate B if not given beforehand
      [B,s,v]=svd(XXt,0);
      B=B(:,1:S);
   end
   if all(size(C)~=[size(X,3) S]) % Generate C if not given beforehand
      % Fixed elements in C
      C=ones(K,S)+randn(K,S)/10;
      if length(Constraints)>5
         C(find(~Cindi))=0;
      end
   end
   if all(size(H)~=[R S]) % Generate H if not given beforehand
      H = rand(R,S);
   end
elseif initi==2
   if Options(5)==0
      disp(' Random initialization')
   end
   if all(size(A)~=[size(X,1) R]) % Generate A if not given beforehand
      A = rand(size(X,1),R);
   end
   if all(size(B)~=[size(X,2) S]) % Generate B if not given beforehand
      B = rand(size(X,2),S);
   end
   if all(size(C)~=[size(X,3) S]) % Generate C if not given beforehand
      C = rand(K,S);
      % Fixed elements in C
      C=ones(K,S)+randn(K,S)/10;
      if length(Constraints)>5
         C(find(~Cindi))=0;
      end
   end
   if all(size(H)~=[R S]) % Generate H if not given beforehand
      H = rand(R,S);
   end
elseif initi==3   
   if Options(5)==0
      disp(' Using input values for initialization')
   end
else
   error(' Options(2) wrongly specified')
end

if initi~=1
   XtX=X(:,:,1)*X(:,:,1)'; % Calculate for evaluating fit (but if initi = 1 it has been calculated)
   for k = 2:K
      XtX = XtX + X(:,:,k)*X(:,:,k)';
   end
end  
fit      = sum(diag(XtX));
olderfit = fit*3;
oldfit   = fit*2;
fit0     = fit;
it       = 0;
Delta    = 1;


if Options(5)==0
   disp(' ')
   disp(' Fitting model ...')
   disp(' Loss-value      Iteration     %VariationExpl')
end

if length(Constraints)>4 %fixed elements
   C=C.*Cindi;
end

% Special feature introduced for imposing specific equality constraints on C in the original application of PARALIND if length(Options)>9
%disp('Set Options(10)=100 for a start')

%lambda = Options(10);if R == 3,   Meq = [1 1 -1 -1 0 0;0 0 1 1 -1 -1];elseif R==2,   Meq = [1 1 -1 -1];,else,   error(' Options not set correct (it''s too long)'),end,MtM = Meq'*Meq*lambda^2;


% Iterative part
while abs(fit-olderfit)>olderfit*ConvCrit & it<MaxIt & fit>1000*eps
   
   olderfit = oldfit;
   oldfit   = fit;
   it       = it + 1;
   
   % Do linesearch
   if rem(it,5)==0&it>10
      %�ndre til 30
      [A,B,C,H,delta]=ParalinLinesearch(A,H,B,C,Ao,Ho,Bo,Co,X,K,MissingElements,MissingOnes,fit,delta);
   end
   Ao=A;Bo=B;Ho=H;Co=C;
   
   
   % Update H
   if Constraints(2)~=-1 % H not fixed
      ZtZ = kron((B'*B).*(C'*C),A'*A);
      XtZ = A'*X(:,:,1)*B*diag(C(1,:));
      for k = 2:K
         XtZ = XtZ + A'*X(:,:,k)*B*diag(C(k,:));
      end       
      if Constraints(2)==0
         h = pinv(ZtZ)*XtZ(:);
      elseif Constraints(2)==1
         h = fastnnls(ZtZ,XtZ(:));
      end
      H = reshape(h,R,S);
   end
   if dbg
      [fit,X] = pflfit(X,A,H,B,C,K,MissingElements,MissingOnes);
      disp(['H:',num2str(fit)]);
   end
   
   % Update A
   Aold = A;
   if Constraints(1)~=-1 % A not fixed
      XtZ = X(:,:,1)*B*diag(C(1,:))*H';
      for k = 2:K
         XtZ = XtZ + X(:,:,k)*B*diag(C(k,:))*H';
      end
      if Constraints(1)==0
         A = XtZ*pinv( H*((B'*B).*(C'*C))*H' );
      elseif Constraints(1)==3 %by Lu, put orthogolity
         A = XtZ*pinv(sqrtm(XtZ'*XtZ)); %by Lu, put orthogolity
      elseif Constraints(1)==1
         ZtZ = H*((B'*B).*(C'*C))*H';
         for i=1:size(A,1)
            A(i,:) = fastnnls(ZtZ,XtZ(i,:)')';
         end
      elseif Constraints(1)==2
         ZtZ = H*((B'*B).*(C'*C))*H';
         A = unimodal(ZtZ,XtZ',A);
      end
   end
   %Check no rank-reduction
   A = chkrank(A,Aold,dbg);
   if dbg
      [fit,X] = pflfit(X,A,H,B,C,K,MissingElements,MissingOnes);
      disp(['A:',num2str(fit)]);
   end
   
   
   % Update B
   Bold = B;
   if Constraints(3)~=-1 % B not fixed
      XtZ = X(:,:,1)'*A*H*diag(C(1,:));
      for k = 2:K
         XtZ = XtZ + X(:,:,k)'*A*H*diag(C(k,:));
      end
      if Constraints(3)==0
         B = XtZ*pinv( (H'*A'*A*H).*(C'*C) );
      elseif Constraints(3)==3  %by Lu, put orthogonality
         B = XtZ*pinv(sqrtm(XtZ'*XtZ)); %by Lu, put orthogolity
      elseif Constraints(3)==1
         ZtZ = (H'*A'*A*H).*(C'*C);
         for i=1:size(B,1)
            B(i,:) = fastnnls(ZtZ,XtZ(i,:)')';
         end
      elseif Constraints(3)==2
         ZtZ = (H'*A'*A*H).*(C'*C);
         B = unimodal(ZtZ,XtZ',B);
      end
      B = chkrank(B,Bold,dbg);   %Check no rank-reduction
      for s = 1:S  % Normalize
         B(:,s) = B(:,s)/norm(B(:,s));
      end
   end
   if dbg
      [fit,X] = pflfit(X,A,H,B,C,K,MissingElements,MissingOnes);
      disp(['B:',num2str(fit)]);
   end

   % Update C
   Cold = C;
   if length(Options)>9 % Just an add-on for imposing specific equality constraint for FIA data
      mtm = MtM*norm(X(:))/(norm(A)*norm(B)*norm(H));
      ZtZ = (B'*B).*(H'*A'*A*H)+mtm;
   else
      ZtZ = (B'*B).*(H'*A'*A*H);
   end
   
   if Constraints(4)~=-1 % C not fixed
      if Constraints(4)==0
         if length(Constraints)==4 % if length(Constraints)>4 fixed to zero on some C parameters
            pZtZ = pinv(ZtZ);
            for k = 1:K
               C(k,:) = (pZtZ*diag(H'*A'*X(:,:,k)*B))';
            end
         else
            for k = 1:K
               id = find(Cindi(k,:));
               XtZ = diag(H'*A'*X(:,:,k)*B);
               C(k,id) = (pinv(ZtZ(id,id))*XtZ(id))';
            end
         end
         
      elseif Constraints(4)==1
         if length(Constraints)>4
            error('Notimplemented')
         end
         for k=1:K
            C(k,:) = fastnnls(ZtZ,diag(H'*A'*X(:,:,k)*B))';
         end
         
      elseif Constraints(4)==2
         %          if length(Constraints)<5
         %             XtZ =[];
         %             for k = 1:K
         %                XtZ = [XtZ;diag(H'*A'*X(:,:,k)*B)'];
         %             end
         %             C = unimodal(ZtZ,XtZ',C);
         %          else
         ZtZ = (B'*B).*(H'*A'*A*H);
         for s = 1:6
            sumC=mean(C')'*2;
            id = find(Cindi(:,s));
            XtZ =[];
            for k = 1:K
               Xk = X(:,:,k)-A*H(:,[1:s-1 s+1:6])*diag(C(k,[1:s-1 s+1:6]))*B(:,[1:s-1 s+1:6])';
               XtZ = [XtZ;diag(H'*A'*Xk*B)'];
            end
            c = (pinv(ZtZ(s,s))*XtZ(id,s)')';
            if s==1|s==3|s==5
               goal = sumC-C(:,s+1);
            else
               goal = sumC-C(:,s-1);
            end
            goal = goal(id);
            p=(1/(1 + lambda))*(c + lambda*goal);
            C(id,s) = ulsr(p,1);
         end
         %         end
      end
   end
   
   %Check no rank-reduction
   C = chkrank(C,Cold,dbg);
   
   if dbg
      [fit,X] = pflfit(X,A,H,B,C,K,MissingElements,MissingOnes); 
      disp(['C:',num2str(fit)]);
   end
   
   % Calculate fit
   [fit,X] = pflfit(X,A,H,B,C,K,MissingElements,MissingOnes);
   
   % Print interim result
   if rem(it,ShowFit)==0|it == 1
      if Options(5)==0
         fprintf(' %12.10f       %g        %3.4f \n',fit,it,100*(1-fit/fit0));
         subplot(2,2,1)
         plot(A),title('First mode')
         subplot(2,2,2)
         plot(B),title('Second mode')
         subplot(2,2,3)
         plot(C),title('Third mode')
         drawnow
      end
   end
   
end

if rem(it,ShowFit)~=0 %Show final fit if not just shown
   if Options(5)==0
      fprintf(' %12.10f       %g        %3.4f \n',fit,it,100*(1-fit/fit0));
   end
end

explainvar=100*(1-fit/fit0);% added by Lu, to show the explained variance



function [fit,X] = pflfit(X,A,H,B,C,K,MissingElements,MissingOnes);

% Calculate fit and impute missing elements from model

fit = 0;
for k = 1:K
   M   = A*H*diag(C(k,:))*B';
   % if missing values replace missing elements with model estimates
   if nargout == 2 
      if any(MissingOnes(:,:,k))
         x=X{k};
         x(find(MissingOnes(:,:,k))) = M(find(MissingOnes(:,:,k)));
         X{k} = x;
      end
   end
   fit = fit + sum(sum(abs (X(:,:,k) - M ).^2));
end


function [b,xi] = fastnnls(x,y,tol,b,eqconst,xi,nnconst)


[m,n] = size(x);
if (nargin < 3 || isempty(tol) || tol == 0 )
  tol = max(size(x))*norm(x,1)*eps;
end
if nargin < 4 || isempty(b);
  b = ones(n,size(y,2));
end
if size(b,2)==1 & size(y,2)>1;  %is b a column vector but y a matrix?
  b(:,2:size(y,2)) = b(:,1);    %copy to every column to match y.
end

if nargin<5 || isempty(eqconst);
  eqconst = zeros(size(b))*nan;  %default is no constraints
end
if size(eqconst,2)==1 & size(b,2)>1;  %is eqconst a column vector but b a matrix?
  eqconst(:,2:size(b,2)) = eqconst(:,1);    %copy to every column to match b.
end
if size(eqconst,1)==1 & size(b,1)>1;  %is eqconst a row vector but b a matrix?
  eqconst(2:size(b,1),:) = eqconst(1,:);    %copy to every row to match b.
end

if nargin<6;  
  xi = [];
  nnconst = [];
else
  if isstruct(xi);
    nnconst = [];
  else
    nnconst = xi;
    xi = [];
  end
end
if isempty(xi); 
    xi = struct('Cache',[]);
end
if isempty(nnconst)
  nnconst = ones(size(b));
  end

if ~isstruct(xi) | checkmlversion('<','6.5') | n>10;
  fn = @proj_old;  %use plain (standard) projection
else
  fn = @fastnnls_proj;  %use caching projection
end
% fn = @proj_robust;

%loop across y columns
y_all = y;
b_all = b;
const_all = eqconst;
for col = 1:size(y,2);
  
  y = y_all(:,col);
  b = b_all(:,col);
  
  eqconst = const_all(:,col);
  noneqconstmap = isnan(eqconst);   %map where 1=non-equality constrained factor
  noneqconstind = find(noneqconstmap);  %lookup table of non-equality constrained factors
  if any(~noneqconstmap);
    y = y - x(:,~noneqconstmap)*eqconst(~noneqconstmap);
    b(~noneqconstmap) = 0;  %set their weight to be zero
  end
  
  p    = (noneqconstmap & b>0)';    %variables which are NOT held at zero
  r    = ~p;                    %variables which ARE held at zero
  b(r) = 0;
  
  [sp,xi] = feval(fn,x,xi,p,y);   %do one projection
  b(p) = sp;   %select reg coef for those factors which were not controlled
  while min(sp) < 0
    b(b<0) = 0;  %assign a zero
    p = (noneqconstmap & b>0)';    %redetermine controlled and uncontrolled vars
    r = ~p;
    [sp,xi] = feval(fn,x,xi,p,y);
    b(p) = sp;
  end
  
  w = x'*(y-x*b);   %correlation beteween x and residuals
  [wmax,ind] = max(w(noneqconstmap));
  ind = noneqconstind(ind);  %locate actual index in unconstrained index list
  flag = 0;
  inloop = 0;
  while (wmax > tol & any(r(noneqconstmap)))
    p(ind) = 1;     %allow that given index to be free
    r(ind) = 0;
    [sp,xi] = feval(fn,x,xi,p,y);
    while (min(sp) < 0) & any(p)   %while any are negative and uncontrolled
      tsp    = zeros(n,1);
      tsp(p) = sp;  
      fb     = (b~=0);
      nrm    = (b(fb)-tsp(fb));
      nrm(nrm<0) = inf;
      rat    = b(fb)./nrm;
      alpha  = min(rat(rat>0));
      alpha  = min([alpha 1]);      %limit to 1
      b = b + alpha*(tsp-b);
      p = (b > tol)';
      r = ~p;
      b(r) = 0;
      [sp,xi] = feval(fn,x,xi,p,y);
    end
    b(p) = sp;
    w = x'*(y-x*b);   %correlation beteween x and residuals
    [wmax,ind] = max(w(noneqconstmap));
    ind = noneqconstind(ind);  %locate actual index in unconstrained index list
    inloop = inloop+1;
    if inloop>100;
      break
    end
    if p(ind)    %already free or stuck in iterations? just leave
      wmax = 0;
    end
  end
  
  if any(~noneqconstmap);
    b(~noneqconstmap) = eqconst(~noneqconstmap);
  end
  
  b_all(:,col) = b;  %store this column's result

  %every n'th point, do a drawnow to allow control-c
  switch mod(col,200)
    case 0
      drawnow;
  end
  
end

b = b_all;

%----------------------------
function [sp,xi] = proj_old(x,xi,p,y);
sp = x(:,p)\y;


%----------------------------
function [sp,xi] = proj_robust(x,xi,p,y);
% sp = x(:,p)\y;
if any(p)
  result = ltsregres(x(:,p),y,'plots',0,'intercept',0);
  sp = result.slope;
else
  sp = [];
end






function B=unimodal(XtX,XtY,Bold)

% Solves the problem min|Y-XB'| subject to the columns of 
% B are unimodal and nonnegative. The algorithm is iterative and
% only one iteration is given, hence the solution is only improving 
% the current estimate
%
% Copyright 1997
%
% Rasmus Bro
% Royal Veterinary & Agricultural University
% Denmark
% rb@kvl.dk
%
% Reference
% Bro and Sidiropoulos, "Journal of Chemometrics", 1998, 12, 223-247. 


B=Bold;
F=size(B,2);
for f=1:F
   xty = XtY(f,:)-XtX(f,[1:f-1 f+1:F])*B(:,[1:f-1 f+1:F])';
   beta=pinv(XtX(f,f))*xty;
   B(:,f)=ulsr(beta',1);
end


function [b,All,MaxML]=ulsr(x,NonNeg);

% ------INPUT------
%
% x          is the vector to be approximated
% NonNeg     If NonNeg is one, nonnegativity is imposed
%
%
%
% ------OUTPUT-----
%
% b 	     is the best ULSR vector
% All 	     is containing in its i'th column the ULSRFIX solution for mode
% 	     location at the i'th element. The ULSR solution given in All
%            is found disregarding the i'th element and hence NOT optimal
% MaxML      is the optimal (leftmost) mode location (i.e. position of maximum)
%
% ___________________________________________________________
%
%
%               Copyright 1997
%
% Nikos Sidiroupolos
% University of Maryland
% Maryland, US
%
%       &
%
% Rasmus Bro
% Royal Veterinary & Agricultural University
% Denmark
%
% 
% ___________________________________________________________


% This file uses MONREG.M

x=x(:);
I=length(x);
xmin=min(x);
if xmin<0
   x=x-xmin;
end


% THE SUBSEQUENT 
% CALCULATES BEST BY TWO MONOTONIC REGRESSIONS

% B1(1:i,i) contains the monontonic increasing regr. on x(1:i)
[b1,out,B1]=monreg(x);

% BI is the opposite of B1. Hence BI(i:I,i) holds the monotonic
% decreasing regression on x(i:I)
[bI,out,BI]=monreg(flipud(x));
BI=flipud(fliplr(BI));

% Together B1 and BI can be concatenated to give the solution to
% problem ULSR for any modloc position AS long as we do not pay
% attention to the element of x at this position


All=zeros(I,I+2);
All(1:I,3:I+2)=B1;
All(1:I,1:I)=All(1:I,1:I)+BI;
All=All(:,2:I+1);
Allmin=All;
Allmax=All;
% All(:,i) holds the ULSR solution for modloc = i, disregarding x(i),


iii=find(x>=max(All)');
b=All(:,iii(1));
b(iii(1))=x(iii(1));
Bestfit=sum((b-x).^2);
MaxML=iii(1);
for ii=2:length(iii)
   this=All(:,iii(ii));
   this(iii(ii))=x(iii(ii));
   thisfit=sum((this-x).^2);
   if thisfit<Bestfit
      b=this;
      Bestfit=thisfit;
      MaxML=iii(ii);
   end
end

if xmin<0
   b=b+xmin;
end


% Impose nonnegativity
if NonNeg==1
   if any(b<0)
      id=find(b<0);
      % Note that changing the negative values to zero does not affect the
      % solution with respect to nonnegative parameters and position of the
      % maximum.
      b(id)=zeros(size(id))+0;
   end
end

function [b,B,AllBs]=monreg(x);

% Monotonic regression according
% to J. B. Kruskal 64
%
% b     = min|x-b| subject to monotonic increase
% B     = b, but condensed
% AllBs = All monotonic regressions, i.e. AllBs(1:i,i) is the 
%         monotonic regression of x(1:i)
%
%
% Copyright 1997
%
% Rasmus Bro
% Royal Veterinary & Agricultural University
% Denmark
% rb@kvl.dk
%


I=length(x);
if size(x,2)==2
   B=x;
else
   B=[x(:) ones(I,1)];
end

AllBs=zeros(I,I);
AllBs(1,1)=x(1);
i=1;
while i<size(B,1)
   if B(i,1)>B(min(I,i+1),1)
      summ=B(i,2)+B(i+1,2);
      B=[B(1:i-1,:);[(B(i,1)*B(i,2)+B(i+1,1)*B(i+1,2))/(summ) summ];B(i+2:size(B,1),:)];
      OK=1;
      while OK
         if B(i,1)<B(max(1,i-1),1)
            summ=B(i,2)+B(i-1,2);
            B=[B(1:i-2,:);[(B(i,1)*B(i,2)+B(i-1,1)*B(i-1,2))/(summ) summ];B(i+1:size(B,1),:)];
            i=max(1,i-1);
         else
            OK=0;
         end
      end
      bInterim=[];
      for i2=1:i
         bInterim=[bInterim;zeros(B(i2,2),1)+B(i2,1)];
      end
      No=sum(B(1:i,2));
      AllBs(1:No,No)=bInterim;
   else
      i=i+1;
      bInterim=[];
      for i2=1:i
         bInterim=[bInterim;zeros(B(i2,2),1)+B(i2,1)];
      end
      No=sum(B(1:i,2));
      AllBs(1:No,No)=bInterim;
   end
end

b=[];
for i=1:size(B,1)
   b=[b;zeros(B(i,2),1)+B(i,1)];
end


function A = chkrank(A,Aold,dbg);

if any(abs(sum(A))<1000*eps)
   A = A*.1+Aold*.9;
   if dbg
      disp(' Rank problems')
   end   
end

function [A,B,C,H,delta]=ParalinLinesearch(A,H,B,C,Ao,Ho,Bo,Co,X,K,MissingElements,MissingOnes,fit,delta);

stop=0;

[fitnew,X] = pflfit(X,A+delta*(A-Ao),H+delta*(H-Ho),B+delta*(B-Bo),C+delta*(C-Co),K,MissingElements,MissingOnes);
while ~stop
   if fitnew>fit
      delta = delta*.6;
      [fitnew,X] = pflfit(X,A+delta*(A-Ao),H+delta*(H-Ho),B+delta*(B-Bo),C+delta*(C-Co),K,MissingElements,MissingOnes);
   end
   if delta<1
      stop=1;
      delta=1;
   end
   if fitnew<fit
      fitnewold=fitnew;
      deltaold=delta;
      delta=delta*1.2;
      [fitnew,X] = pflfit(X,A+delta*(A-Ao),H+delta*(H-Ho),B+delta*(B-Bo),C+delta*(C-Co),K,MissingElements,MissingOnes);
      if fitnew>fit|fitnew>fitnewold
         stop = 1;
         delta=deltaold;
         A=A+delta*(A-Ao);
         H=H+delta*(H-Ho);
         B=B+delta*(B-Bo);
         C=C+delta*(C-Co);
      end
   end
end