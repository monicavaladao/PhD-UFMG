function [XMIN,OBJMIN,G,GXMIN,GOBJMIN,GNfobj,numfobj,RUNTime] = algoritmo_de_best_normalizado(fobj,P,jP,lb,ub,N,countnumfobj,numfobj)
% ALGORITMO_DE_NORMALIZADO: Implementa o algoritmo DE/current-to-best/1
% malizado no intervalo [0,1]
% Entrada:
%   fobj       : Fun��o objetivo
%   P,jP       : Popula��o corrente e respectiva avalia��o
%   [lb,ub]    : Limites para as vari�veis de decis�o
%   n,N        : N�mero de vari�veis e tamanho de P
%   countnumfobj: N�mero de avalia��es j� realizadas
%   numfobj    : N�mero m�ximo de avalia��es em fobj
% Saida:
%  P,jP        : Popula��o e respectiva avalia��o atualizada
%  XMIN,OBJMIN : Melhor solu��o e respectiva avalia��o
%  G           : N�mero de gera��e realizadas
% GXMIN,GOBJMIN: Melhor solu��o e respectiva avalia��o a cada itera��o
% GNfobj       : N�mero de avalia��es realizadas a cada itera��o
%  numobj      : N�mero de avalia��es realizadas

% Parameters
g = 1;          % Generation counter

C  = 0.5;       % Crossover rate

t0_start = cputime;

t0_end = 0;

% Stop criterion
while (countnumfobj < numfobj)
    
    % Store the best solution from the current generation
    [GOBJMIN(g),j_best] = min(jP);
    GXMIN(:,g) = (ub(:) - lb(:)).*P(:,j_best) + lb(:);% Solu��o desnormalizada
    G(g) = g;
    GNfobj(g) = countnumfobj;
	
	RUNTime(g) = t0_end;
    
    % Differential mutataion
    V = mutation(P,j_best);
        
    % Discrete recombination
    U = recombination(P,V,C);
    
    % Truncation and Evaluation
    U  = truncation(U);
    
    % Avalia as novas solu��es em P
    if countnumfobj + N >= numfobj
        % Identifica o n�mero de avalia��es que ainda restam
        lastsol = numfobj - countnumfobj;
        [U,jU] = seleciona_avalia_solucoes_aleatorias(fobj,P,jP,U,N,lb,ub,lastsol);
        % Atualiza contador de avalia��es
        countnumfobj = countnumfobj + lastsol;
    else
        % Avalia as N novas solu��es
        jU = evaluation(fobj,U,lb,ub);
        % Atualiza contador de avalia��es
        countnumfobj = countnumfobj + N;
    end
    
    % Survival selection
    [P,jP] = selection(P,jP,U,jU);
    
    % Plot the best current solution
    %plotter(fobj,lb,ub,g,P,jP);
	
	% Time
	t0_end = cputime - t0_start;
       
    % Update counter
    g = g + 1;
    
end

% Store the best solution from the current generation
[GOBJMIN(g),j_best] = min(jP);
GXMIN(:,g) = (ub(:) - lb(:)).*P(:,j_best) + lb(:);% Solu��o desnormalizada
G(g) = g;
GNfobj(g) = countnumfobj;

RUNTime(g) = t0_end;

% Retorma a melhor solu��o
[OBJMIN,j_best] = min(jP);
XMIN = P(:,j_best);
% Desnormaliza a solu��o
XMIN = (ub(:) - lb(:)).*XMIN + lb(:);

% figure
% plot(b,'k-o','LineWidth',2,'MarkerSize',2)
% title('Differential Evolution Algorithm')
% xlabel('generation')
% ylabel('f(x) (best per generation)')

end

% Fun��es auxiliares
%(1) % Fun��o de muta��o
function V = mutation(P,jbest)

%
% DE/rand/1/bin
%

[n,N] = size(P);
V     = zeros(n,N);

for i = 1:N
    
    j = randperm(N);
    while (sum(i == j(1:2)) > 0)
        j = randperm(N);
    end
    
    % Random weight
    F = 0.4 + 0.5*rand();
    V(:,i) = P(:,jbest) + F*(P(:,j(1))-P(:,j(2)));
end
end
%
%(2) % Fun��o de recombina��o
function U = recombination(P,V,C)

[n,N] = size(P);
U = zeros(n,N);

for i = 1:N
    
    deltai = round((n-1)*rand(1)+1);
    
    for j = 1:n
        if (rand <= C || j == deltai)
            U(j,i) = V(j,i);            
        else            
            U(j,i) = P(j,i);
        end
    end
end
end