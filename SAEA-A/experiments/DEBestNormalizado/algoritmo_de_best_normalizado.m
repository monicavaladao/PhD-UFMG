function [XMIN,OBJMIN,G,GXMIN,GOBJMIN,GNfobj,numfobj,RUNTime] = algoritmo_de_best_normalizado(fobj,P,jP,lb,ub,N,countnumfobj,numfobj)
% ALGORITMO_DE_NORMALIZADO: Implementa o algoritmo DE/current-to-best/1
% malizado no intervalo [0,1]
% Entrada:
%   fobj       : Função objetivo
%   P,jP       : População corrente e respectiva avaliação
%   [lb,ub]    : Limites para as variáveis de decisão
%   n,N        : Número de variáveis e tamanho de P
%   countnumfobj: Número de avaliações já realizadas
%   numfobj    : Número máximo de avaliações em fobj
% Saida:
%  P,jP        : População e respectiva avaliação atualizada
%  XMIN,OBJMIN : Melhor solução e respectiva avaliação
%  G           : Número de geraçõe realizadas
% GXMIN,GOBJMIN: Melhor solução e respectiva avaliação a cada iteração
% GNfobj       : Número de avaliações realizadas a cada iteração
%  numobj      : Número de avaliações realizadas

% Parameters
g = 1;          % Generation counter

C  = 0.5;       % Crossover rate

t0_start = cputime;

t0_end = 0;

% Stop criterion
while (countnumfobj < numfobj)
    
    % Store the best solution from the current generation
    [GOBJMIN(g),j_best] = min(jP);
    GXMIN(:,g) = (ub(:) - lb(:)).*P(:,j_best) + lb(:);% Solução desnormalizada
    G(g) = g;
    GNfobj(g) = countnumfobj;
	
	RUNTime(g) = t0_end;
    
    % Differential mutataion
    V = mutation(P,j_best);
        
    % Discrete recombination
    U = recombination(P,V,C);
    
    % Truncation and Evaluation
    U  = truncation(U);
    
    % Avalia as novas soluções em P
    if countnumfobj + N >= numfobj
        % Identifica o número de avaliações que ainda restam
        lastsol = numfobj - countnumfobj;
        [U,jU] = seleciona_avalia_solucoes_aleatorias(fobj,P,jP,U,N,lb,ub,lastsol);
        % Atualiza contador de avaliações
        countnumfobj = countnumfobj + lastsol;
    else
        % Avalia as N novas soluções
        jU = evaluation(fobj,U,lb,ub);
        % Atualiza contador de avaliações
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
GXMIN(:,g) = (ub(:) - lb(:)).*P(:,j_best) + lb(:);% Solução desnormalizada
G(g) = g;
GNfobj(g) = countnumfobj;

RUNTime(g) = t0_end;

% Retorma a melhor solução
[OBJMIN,j_best] = min(jP);
XMIN = P(:,j_best);
% Desnormaliza a solução
XMIN = (ub(:) - lb(:)).*XMIN + lb(:);

% figure
% plot(b,'k-o','LineWidth',2,'MarkerSize',2)
% title('Differential Evolution Algorithm')
% xlabel('generation')
% ylabel('f(x) (best per generation)')

end

% Funções auxiliares
%(1) % Função de mutação
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
%(2) % Função de recombinação
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