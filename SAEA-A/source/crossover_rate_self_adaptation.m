% function [CR, CRm, S_A] = crossover_rate_self_adaptation(CR, CRm, S_A, N, n, g)
% % CROSSOVER_RATE_ADAPTED: Perform a self-adaptation on crossover rate.
% % Input:
% %   CR: Vector crossover rate (Nx1)
% %   CRm: Parameter to update each CRi
% %   S_A: Structure whith informations at a given generation
% %       CRrec: CR values associated with offspring successfully entering in
% %       the each generation
% %       Pimp: Vector whith the diferences on percentage of improvement
% %       related to CRrec
% %   N: Population size
% %   n: Number of variables
% %   g: Current counter of the number times in which the DE operator is used
% % Output:
% %   CR: Vector crossover rate (Nx1)
% %   CRm: Parameter to update each CRi
% %   S_A: Structure whith informations at a given generation
% %       CRrec: CR values associated with offspring successfully entering in
% %       the each generation
% %       Pimp: Vector whith the diferences on percentage of improvement
% %       related to CRrec
% 
% % Number of solutions corresponding to vector CR
% %N = length(CR);
% 
% % Define parameters according the number of variables n
% if n <= 100
%     Ng = 5;
%     Nmod = 5;
% else
%     Ng = 10;
%     Nmod = 10;
% end
% 
% %if g <= 5
% if g <= Ng
%     % Keep the initial vector crossover rate
%     CR = CR;
% else
%     % Sef-adaption on crossover rate
%     if mod(g,20) == 0
%     %if mod(g,10) == 0
%         Rec_CR = [S_A.CRrec];
%         Rec_Pimp = [S_A.Pimp];
%         if isempty(Rec_CR) == 0 && isempty(Rec_Pimp) == 0
%             W = Rec_Pimp./sum(Rec_Pimp);
%             CRm = sum(W.*Rec_CR);
%             
%             % Upadate CR
%             if isnan(CRm)
%                 CR = CR;
%             elseif isinf(CRm)
%                 CR = CR;
%             else
%                 [CR] = crossover_rate(CRm, N);
%             end
%             
%             % Update S_A structure
%             S_A.CRrec = [];
%             S_A.Pimp = [];
%         end
%     else
%         %if mod(g,5) == 0
%         if mod(g,Nmod) == 0
%             % Upadate CR
%             [CR] = crossover_rate(CRm, N);
%         end
%     end
% 
% end
% 
% 
% end

% % MELHOR CROSSOVER UPDTATE ATÉ AGORA
% function [CR, CRm, S_A] = crossover_rate_self_adaptation(CR, CRm, S_A, N, n, g)
% % CROSSOVER_RATE_ADAPTED: Perform a self-adaptation on crossover rate.
% % Input:
% %   CR: Vector crossover rate (Nx1)
% %   CRm: Parameter to update each CRi
% %   S_A: Structure whith informations at a given generation
% %       CRrec: CR values associated with offspring successfully entering in
% %       the each generation
% %       Pimp: Vector whith the diferences on percentage of improvement
% %       related to CRrec
% %   N: Population size
% %   n: Number of variables
% %   g: Current counter of the number times in which the DE operator is used
% % Output:
% %   CR: Vector crossover rate (Nx1)
% %   CRm: Parameter to update each CRi
% %   S_A: Structure whith informations at a given generation
% %       CRrec: CR values associated with offspring successfully entering in
% %       the each generation
% %       Pimp: Vector whith the diferences on percentage of improvement
% %       related to CRrec
% 
% % Number of solutions corresponding to vector CR
% %N = length(CR);
% 
% % Define parameters according the number of variables n
% if n <= 100
%     Ng = 5;
%     Nmod = 5;
% else
%     Ng = 10;
%     Nmod = 10;
% end
% 
% aux_CRm = CRm;%0.5;
% 
% %if g <= 5
% if g <= Ng
%     % Keep the initial vector crossover rate
%     CR = CR;
% else
%     % Sef-adaption on crossover rate
%     if mod(g,20) == 0
%         Rec_CR = [S_A.CRrec];
%         Rec_Pimp = [S_A.Pimp];
%         if isempty(Rec_CR) == 0 && isempty(Rec_Pimp) == 0
%             W = Rec_Pimp./sum(Rec_Pimp);
%             CRm = sum(W.*Rec_CR);
%             % Upadate CR
%             if isnan(CRm)
%                 CRm = aux_CRm;
%             elseif isinf(CRm)
%                 CRm = aux_CRm;
%             elseif CRm == 0
%                 CRm = aux_CRm;
%             end
%             [CR] = crossover_rate(CRm, N);
%             % Update S_A structure
%             S_A.CRrec = [];
%             S_A.Pimp = [];
% 
%         end
%     else
% %        if mod(g,Nmod) == 0
%             % Upadate CR
%             [CR] = crossover_rate(CRm, N);
% %        end
%     end
% 
% end
% 
% 
% end


% MELHOR CROSSOVER UPDTATE ATÉ AGORA
function [CR, CRm, S_A] = crossover_rate_self_adaptation(CR, CRm, S_A, N, n, g)
% CROSSOVER_RATE_ADAPTED: Perform a self-adaptation on crossover rate.
% Input:
%   CR: Vector crossover rate (Nx1)
%   CRm: Parameter to update each CRi
%   S_A: Structure whith informations at a given generation
%       CRrec: CR values associated with offspring successfully entering in
%       the each generation
%       Pimp: Vector whith the diferences on percentage of improvement
%       related to CRrec
%   N: Population size
%   n: Number of variables
%   g: Current counter of the number times in which the DE operator is used
% Output:
%   CR: Vector crossover rate (Nx1)
%   CRm: Parameter to update each CRi
%   S_A: Structure whith informations at a given generation
%       CRrec: CR values associated with offspring successfully entering in
%       the each generation
%       Pimp: Vector whith the diferences on percentage of improvement
%       related to CRrec

% Number of solutions corresponding to vector CR
%N = length(CR);

% Define parameters according the number of variables n
if n <= 100
    Ng = 5;
    Nmod = 5;
else
    Ng = 10;
    Nmod = 10;
end

aux_CRm = CRm;%0.5;

%if g <= 5
if g <= Ng
    % Keep the initial vector crossover rate
    CR = CR;
else
    % Sef-adaption on crossover rate
    Rec_CR = [S_A.CRrec];
    Rec_Pimp = [S_A.Pimp];
    if isempty(Rec_CR) == 0 && isempty(Rec_Pimp) == 0
        W = Rec_Pimp./sum(Rec_Pimp);
        CRm = sum(W.*Rec_CR);
        % Upadate CR
        if isnan(CRm)
            CRm = aux_CRm;
        elseif isinf(CRm)
            CRm = aux_CRm;
        elseif CRm == 0
            CRm = aux_CRm;
        end
    end
    [CR] = crossover_rate(CRm, N);
    %if mod(g,20) == 0
    if mod(g,40) == 0
        % Update S_A structure
        S_A.CRrec = [];
        S_A.Pimp = [];    
    end

end


end