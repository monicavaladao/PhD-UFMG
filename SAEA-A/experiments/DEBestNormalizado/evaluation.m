% (5) Função de avaliação

function jX = evaluation(fobj,X,lb,ub)

N  = size(X,2);
jX = nan(1,N);

% Desnormalize population
X = (repmat(ub(:),1,N)-repmat(lb(:),1,N)).*X + repmat(lb(:),1,N);

for i = 1:N
    jX(i) = feval(fobj,X(:,i)'); 
end 
% 
% switch fobj
%     case 'LS_fun'
%         for i = 1:N
%             jX(i) = LS_fun(X(:,i)');
%         end
%     otherwise
%         for i = 1:N
%             jX(i) = feval(fobj,X(:,i));
%         end  
% end
       

end