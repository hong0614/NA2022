function [a,s]=dls(x,y,n)
[~,col]=size(x);
G=zeros(n);
c=zeros(n,1);

for i=1:n
    for j=1:n
        if(j==1&&i==1)
            G(1,1)=col;
        else
            G(i,j)=sum(x.^(i-1).*x.^(j-1));
        end
    end
end

for i = 1:n
    c(i,1)=sum(x.^(i-1).*y);
end

a=G\c;
s=cond(G,2);

