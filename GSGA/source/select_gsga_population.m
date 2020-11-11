function [pop_X,pop_y] = select_gsga_population(pop_X,pop_y,New_P,New_y,N3)

[~,idx_sort] = sort(pop_y,'ascend');
pop_X = pop_X(idx_sort,:);
pop_y = pop_y(idx_sort);

aux_P = pop_X(1:N3,:);
aux_y = pop_y(1:N3);
pop_X = [];
pop_y = [];
pop_X = [aux_P;New_P];
pop_y = [aux_y(:);New_y(:)];


end