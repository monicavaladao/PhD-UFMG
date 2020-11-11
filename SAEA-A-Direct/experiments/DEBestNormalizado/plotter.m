
function plotter(fobj,nomealgoritmo,lb,ub,g,P,jP)

[s,N] = size(P);

% Desnormalize population
P = (repmat(ub(:),1,N)-repmat(lb(:),1,N)).*P + repmat(lb(:),1,N);

if (s == 2)
    [x,y] = meshgrid (lb(1):1/4:ub(1), lb(2):1/4:ub(2));
    z = eval_matrix(fobj,x,y);
    contour(x,y,z); 

    hold on    
    plot(P(1,:), P(2,:), 'k*');
    title(nomealgoritmo)
    xlabel('x_1')
    ylabel('x_2')
    grid on
    hold off
end

% Print the result at each iteration
[opt,iopt]  = min(jP);
fprintf(1, 'f(x) = ')
fprintf(1, '%+6.4f  ', opt)
fprintf(1, '\t x = [')
fprintf(1, '%+6.4f  ', P(:,iopt)')
fprintf(1, '\b\b]; \t g = %d\n', g)

drawnow



function z = eval_matrix(fobj,x,y)

[u,v] = size(x);
z = zeros(u,v);

for i = 1:u
    for j = 1:v
        z(i,j) = feval(fobj,[x(i,j) y(i,j)]);
    end
end

