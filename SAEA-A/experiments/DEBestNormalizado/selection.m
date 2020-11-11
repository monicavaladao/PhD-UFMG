% (4) Função de seleção

function [S,jS] = selection(X,jX,P,jP)

[n,N] = size(P);

S  = zeros(n,N);
jS = zeros(1,N);

i = find(jX<=jP);
j = find(jX>jP);

S(:,i) = X(:,i);
jS(i)  = jX(i);

S(:,j) = P(:,j);
jS(j)  = jP(j);

end