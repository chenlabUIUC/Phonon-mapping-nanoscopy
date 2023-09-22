function d=LatticeError(X,x)
%x=[80,0,0,80,-4,4,-4,4];
%x=[1 ,2,3, 4, 5,6, 7,8,];
d=0;
X2=X;
l=length(X);
if length(X) == length(x(5):x(6))*length(x(7):x(8))
for n=x(5):x(6)
    for m=x(7):x(8)
        [M,I]=min((X2(:,1)-(n.*x(1)+m.*x(3))).^2+(X2(:,2)-(n.*x(2)+m.*x(4))).^2);
        d=d+M/l;
        X2(I,1)=10^10;
        X2(I,2)=10^10;
    end
end
else
    d=10^10;
end