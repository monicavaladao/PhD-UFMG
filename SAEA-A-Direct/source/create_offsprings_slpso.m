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
% SL-PSO version adpated to use in the proposed SAEA

function [P,V,CR] = create_offsprings_slpso(pop_X, pop_y, pop_V, pop_CR, lb, ub, bestever)


%d: Size and dimensionality of pop_X
[m,d] = size(pop_X);

% lu: define the upper and lower bounds of the variables
lu = [lb;ub];

% Identify informations
p = pop_X;
fitness = pop_y;
%v = zeros(m,d);
v = pop_V;
CR = pop_CR;

%parameter initiliaztion
M = m - floor(d/10);
c3 = (d/M)*0.01;
PL = zeros(m,1);

%alf = 0.5;
for i = 1 : m
    PL(i) = (1 - (i - 1)/m)^log(sqrt(ceil(d/M)));
    %PL(i) = (1 - (i - 1)/m)^(alf*log(ceil(d/M)));
end

% population sorting
[fitness rank] = sort(fitness, 'descend');
p = p(rank,:);
v = v(rank,:);
PL = PL(rank);
besty = fitness(m);
bestp = p(m, :);
bestever = min(besty, bestever);

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

% Return offsprings
P = p;
V = v;
% population sorting
[fitness rank] = sort(fitness, 'ascend');
P = P(rank,:);
V = V(rank,:);
CR = CR(rank);
end
