function z=TC(X)%X is a cell, TC means Tucker Congurancy
m=length(X);
for i=1:m(1)
    s=size(X{i});
    for j=1:s(2)-1
        for k=(j+1):s(2)
            y(i,(j-1)*s(2)-(j+1)/2*j+k)=(X{i}(:,j)'*X{i}(:,k))/(norm(X{i}(:,j))*norm(X{i}(:,k)));
        end
    end
end
z=prod(y);


