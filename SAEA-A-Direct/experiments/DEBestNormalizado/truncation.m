% (3) Fun��o para corrigir os limites de vari�veis
function M = truncation(M)

M(M<0) = 0;
M(M>1) = 1;

end