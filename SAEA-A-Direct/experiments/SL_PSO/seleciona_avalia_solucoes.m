function [PLast,FPLast,VLast] = seleciona_avalia_solucoes(fobj,P,FP,V,Pop,PopV,lastsol)
% SELECIONA_AVALIA_SOLUCOES_ALEATORIAS: Função que seleciona lastsol
% soluções aleatórias em Pop e avalia em fobj.
% Entrada:
%   fobj       : Função objetivo
%   P,FP       : População corrente e respectiva avaliação
%   Pop        : População gerada de acordo com o algoritmo escolhido
%   lastsol    : Número permitido de avaliação em fobj
% Saída:
%  PLast,FPLast: População e respectiva avaliação considerando as
%                avaliações de lastgen soluções

% Inicializa matrizes de armazenamento
PLast = P;
FPLast = FP;
VLast = V;

% Escolhe lastsol soluções
Last = 1:lastsol;

% Identifica e avalia as lastgen soluções aleatórias
PopLast = Pop(Last,:);
FPopLast = feval_all(fobj,PopLast);
PopVLast = PopV(Last,:);

% Retorna população e respectiva avaliação adaptada com as lastgen soluções aleatórias
PLast(Last,:) = PopLast;
FPLast(Last) = FPopLast;
VLast(Last,:) = PopVLast;

end