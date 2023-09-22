function x=Evolutionaryoptomize(f,lb,ub,x0,Pop,survive,child,fresh,mutation,Gmax);
q=length(lb)
%generation 0
Population=rand(Pop-1,q).*(ub-lb)-lb;
Population=vertcat(x0,Population);
Eval=zeros(Pop,1);
for i=1:Pop
    Eval(i)=f(Population(i,:));
end
M=min(Eval)
figure(); hold on; 
plot(0,M,'.k'); drawnow;
set(gca, 'YScale', 'log'); %set(gca, 'XScale', 'log');
for G=1:Gmax
    [B,I]=mink(Eval,survive);
    x=Population(I(1),:);
    Population2=x;
    for i=2:survive
        Population2=vertcat(Population2,Population(I(i),:));
    end
    Population2=vertcat(Population2,rand(fresh,q).*(ub-lb)-lb);
    for i=1:child
        Q=rand(1,q);
        k= randi([1 survive+fresh],1,2);
        P=(Population2(k(1),:)+randn(1,q).*mutation).*Q+(Population2(k(2),:)+randn(1,q).*mutation).*(1.-Q);
        Population2=vertcat(Population2,P);
    end
    Population(:,:)=Population2(:,:);
    Eval=zeros(Pop,1);
    for i=1:Pop
        Eval(i)=f(Population(i,:));
    end
    plot(G,min(Eval),'.k'); drawnow;
end