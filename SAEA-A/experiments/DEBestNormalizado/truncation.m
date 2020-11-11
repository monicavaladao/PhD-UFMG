% (3) Função para corrigir os limites de variáveis
function M = truncation(M)

M(M<0) = 0;
M(M>1) = 1;

end