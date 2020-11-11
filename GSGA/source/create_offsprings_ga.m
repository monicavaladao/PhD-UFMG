function [P_child] = create_offsprings_ga(pop_X,RBF_X,lb,ub,n,N,N1,N2,pcross,pmut)
% CREATE_OFFSPRINGS_GA: Use GA operators to create the child population.
% Input:
%   pop_X: Current population
%   pop_y: Value of each solution in pop_X
%   RBF_X: Predicted population points
%   lb: Lower bound
%   ub: Upper bound
%   n: Number of variable
%   N: Number of solutions in pop_X
%   N1: gsga parameter
%   N2: gsga parameter
% Output:
%   P_child: The population child

% Rank selection
pdf = (N +1 -[1:N]')/(N*(N+1)/2);
pdf = pdf/sum(pdf);
cdf = cumsum(pdf);

% Selection of population points and predicted population points
PS = []; % Population points
PPS = []; % Predicted population points
for i = 1:N1 % E se for vazio ???
    r = rand;
    idx = find(r <= cdf);
    PS = [PS;pop_X(idx,:)];
    PPS = [PPS;RBF_X(idx,:)];
end

% Uniform crossover
FP = []; % Father points
MP = []; % Mother points
PFP = []; % Predicted father points
PMP = []; % Predicted mother points

for i = 1:N1
    FP = [FP;PS(i,:)];
    MP = [MP;PS(i+1,:)];
    PFP = [PFP;PPS(i,:)];
    PMP = [PMP;PPS(i+1,:)];
end

% Child population after crossover
Child = [];
for i = 1:(N1/2)
    for j = 1:n
        pc = rand;
        if pc <= pcross
            Child1(i,j) = PFP(i,j);
            Child2(i,j) = MP(i,j);
        else
            Child1(i,j) = MP(i,j);
            Child2(i,j) = PFP(i,j);
        end
    end
end
Child = [Child;Child1;Child2];

% Mutation
for i = 1:N2
    j = randi([1 N1]);    
    for k = 1:n
        rm = rand;
        if rm <= pmut
            Child(j,k) = lb(k) + (ub(k) - lb(k))*rand;
        end
    end
    
end
P_child = Child;
end