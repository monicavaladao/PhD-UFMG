function [PLast,FPLast] = seleciona_avalia_solucoes_aleatorias(fobj,P,FP,Pop,N,lb,ub,lastsol)
% SELECIONA_AVALIA_SOLUCOES_ALEATORIAS: Fun��o que seleciona lastsol
% solu��es aleat�rias em Pop e avalia em fobj.
% Entrada:
%   fobj       : Fun��o objetivo
%   P,FP       : Popula��o corrente e respectiva avalia��o
%   Pop        : Popula��o gerada de acordo com o algoritmo escolhido
%   N          : Tamanho da popula��o P
%   lastsol    : N�mero permitido de avalia��o em fobj
% Sa�da:
%  PLast,FPLast: Popula��o e respectiva avalia��o considerando as
%                avalia��es de lastgen solu��es

% Inicializa matrizes de armazenamento
PLast = P;
FPLast = FP;

% Escolhe lastsol solu��es aleat�rias
randLast = randi(N,1,lastsol);

% Identifica e avalia as lastgen solu��es aleat�rias
PopLast = Pop(:,randLast);
FPoplast = evaluation(fobj,PopLast,lb,ub);

% Retorna popula��o e respectiva avalia��o adaptada com as lastgen solu��es aleat�rias
PLast(:,randLast) = PopLast;
FPLast(:,randLast) = FPoplast;

end