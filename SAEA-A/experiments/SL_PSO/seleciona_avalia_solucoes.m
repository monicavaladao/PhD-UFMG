function [PLast,FPLast,VLast] = seleciona_avalia_solucoes(fobj,P,FP,V,Pop,PopV,lastsol)
% SELECIONA_AVALIA_SOLUCOES_ALEATORIAS: Fun��o que seleciona lastsol
% solu��es aleat�rias em Pop e avalia em fobj.
% Entrada:
%   fobj       : Fun��o objetivo
%   P,FP       : Popula��o corrente e respectiva avalia��o
%   Pop        : Popula��o gerada de acordo com o algoritmo escolhido
%   lastsol    : N�mero permitido de avalia��o em fobj
% Sa�da:
%  PLast,FPLast: Popula��o e respectiva avalia��o considerando as
%                avalia��es de lastgen solu��es

% Inicializa matrizes de armazenamento
PLast = P;
FPLast = FP;
VLast = V;

% Escolhe lastsol solu��es
Last = 1:lastsol;

% Identifica e avalia as lastgen solu��es aleat�rias
PopLast = Pop(Last,:);
FPopLast = feval_all(fobj,PopLast);
PopVLast = PopV(Last,:);

% Retorna popula��o e respectiva avalia��o adaptada com as lastgen solu��es aleat�rias
PLast(Last,:) = PopLast;
FPLast(Last) = FPopLast;
VLast(Last,:) = PopVLast;

end