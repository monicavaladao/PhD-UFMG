%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Implementation of a social learning PSO (SL-PSO) for scalable optimization
%
%  See the details of SL-PSO in the following paper
%  R. Cheng and Y. Jin, A Social Learning Particle Swarm Optimization Algorithm for Scalable Pptimization,
%  Information Sicences, 2014
%
%  The source code SL-PSO is implemented by Ran Cheng 
%
%  If you have any questions about the code, please contact: 
%  Ran Cheng at r.cheng@surrey.ac.uk 
%  Prof. Yaochu Jin at yaochu.jin@surrey.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Adapted in 05/02/2020
function [XMIN,OBJMIN,G,GXMIN,GOBJMIN,GNfobj,numfobj,RUNTime] = sl_pso(fobj,P,FP,lb,ub,countnumfobj,numfobj)

%maxfe: maximal number of fitness evaluations
maxfe = numfobj;

%d: Size and dimensionality of P
[m,d] = size(P);

% lu: define the upper and lower bounds of the variables
lu = [lb;ub];

% Identify informations
p = P;
fitness = FP;

%parameter initilization
M = m - floor(d/10);
c3 = d/M*0.01;
PL = zeros(m,1);

for i = 1 : m
    PL(i) = (1 - (i - 1)/m)^log(sqrt(ceil(d/M)));
end

% initialization
v = zeros(m,d);
bestever = 1e200;

FES = countnumfobj;
g = 1;
    

t0_start = cputime;
t0_end = 0;

% main loop
while(FES < maxfe)
    
    % population sorting
    [fitness, rank] = sort(fitness, 'descend');
    p = p(rank,:);
    v = v(rank,:);
	PL = PL(rank);
    besty = fitness(m);
    bestp = p(m, :);
    bestever = min(besty, bestever);
    
    % Store the best solution from the current generation
    GOBJMIN(g) = besty;
    GXMIN(g,:) = bestp;
    G(g) = g;
    GNfobj(g) = FES;
	RUNTime(g) = t0_end;
    
    % Store current information
    P = p;
    FP = fitness;
    V = v;

    % center position
    center = ones(m,1)*mean(p);

    %random matrix 
    %rand('seed', sum(100 * clock));
    randco1 = rand(m, d);
    %rand('seed', sum(100 * clock));
    randco2 = rand(m, d);
    %rand('seed', sum(100 * clock));
    randco3 = rand(m, d);
    winidxmask = repmat([1:m]', [1 d]);
    winidx = winidxmask + ceil(rand(m, d).*(m - winidxmask));
    pwin = p;
    for j = 1:d
            pwin(:,j) = p(winidx(:,j),j);
    end

    % social learning
     lpmask = repmat(rand(m,1) < PL, [1 d]);
     lpmask(m,:) = 0;
     v1 =  1*(randco1.*v + randco2.*(pwin - p) + c3*randco3.*(center - p));
     p1 =  p + v1;   


     v = lpmask.*v1 + (~lpmask).*v;         
     p = lpmask.*p1 + (~lpmask).*p;

     % boundary control
    for i = 1:m - 1
        p(i,:) = max(p(i,:), lu(1,:));
        p(i,:) = min(p(i,:), lu(2,:));
    end

    % Avalia as novas soluções em P
    if (FES + (m - 1)) >= numfobj
        % Identifica o número de avaliações que ainda restam
        lastsol = numfobj - FES;
        [PLast,FPLast,VLast] = seleciona_avalia_solucoes(fobj,P(1:(m - 1),:),FP(1:(m - 1)),...
                               V(1:(m - 1),:),p(1:(m - 1),:),v(1:(m - 1),:),lastsol);
        p(1:m - 1,:) = PLast;
        fitness(1:m - 1,:) = FPLast(:);
        v(1:m - 1,:) = VLast;
        % Atualiza contador de avaliações
        FES = FES + lastsol;
    else
        % Avalia as (m - 1) novas soluções
        % fitness evaluation
        aux_fitness = feval_all(fobj,p(1:m - 1,:));
        fitness(1:m - 1,:) = aux_fitness(:); 
        % Atualiza contador de avaliações
        FES = FES + m - 1;
    end
    
%    % Print information
%   fprintf('Best fitness: %e\n', min(bestever,min(fitness))); 

    % Time
	t0_end = cputime - t0_start;
    
    g = g + 1;
end

% Store the best solution from the current generation
[GOBJMIN(g),j_best] = min(fitness);
GXMIN(g,:) = p(j_best,:);
G(g) = g;
GNfobj(g) = FES;
RUNTime(g) = t0_end;

% Retorma a melhor solução
[OBJMIN,j_best] = min(GOBJMIN);
XMIN = GXMIN(j_best,:);

    
end


    

