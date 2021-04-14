
% This function is used to compute the Tucker congruency
% X is a cell, 
% z is a vector; 
% z stores the TC values of the j-th (j=1,2...) and k-th (j+1,...) component

function z=TC(X) 

m=length(X);

for i=1:m     %Loop in each mode
    s=size(X{i});
    for j=1:s(2)-1
        for k=(j+1):s(2)
            y(i,(j-1)*s(2)-(j+1)/2*j+k)=(X{i}(:,j)'*X{i}(:,k))/(norm(X{i}(:,j))*norm(X{i}(:,k)));
        end
    end
end
z=prod(y);



