function [e,y] = feval_all_two_output(fobj, X)
% FEVAL_ALL: Evaluate each column in X on fobj.

% Number of entries to evaluate
[N,~] = size(X);

% Evaluate all entries
y = zeros(N, 1);
e = zeros(N, 1);
for i = 1:N
    [e(i,1),y(i,1)] = feval(fobj, X(i,:));
end

end