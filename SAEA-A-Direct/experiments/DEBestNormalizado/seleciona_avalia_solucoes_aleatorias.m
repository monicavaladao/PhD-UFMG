function [PLast,FPLast] = seleciona_avalia_solucoes_aleatorias(fobj,P,FP,Pop,N,lb,ub,lastsol)
% SELECIONA_AVALIA_SOLUCOES_ALEATORIAS: Função que seleciona lastsol
% soluções aleatórias em Pop e avalia em fobj.
% Entrada:
%   fobj       : Função objetivo
%   P,FP       : População corrente e respectiva avaliação
%   Pop        : População gerada de acordo com o algoritmo escolhido
%   N          : Tamanho da população P
%   lastsol    : Número permitido de avaliação em fobj
% Saída:
%  PLast,FPLast: População e respectiva avaliação considerando as
%                avaliações de lastgen soluções

% Inicializa matrizes de armazenamento
PLast = P;
FPLast = FP;

% Escolhe lastsol soluções aleatórias
randLast = randi(N,1,lastsol);

% Identifica e avalia as lastgen soluções aleatórias
PopLast = Pop(:,randLast);
FPoplast = evaluation(fobj,PopLast,lb,ub);

% Retorna população e respectiva avaliação adaptada com as lastgen soluções aleatórias
PLast(:,randLast) = PopLast;
FPLast(:,randLast) = FPoplast;

end