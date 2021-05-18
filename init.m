%% CapSA (Capuchin Swarm Algorithm)
% Citation details:
% Braik, Malik, Alaa Sheta, and Heba Al-Hiary. "A novel meta-heuristic search algorithm for solving 
% optimization problems: capuchin search algorithm.
% " Neural Computing and Applications (2020): 1-33

% Programmed by Malik Braik & Prof. Alaa Sheta
% Al-Balqa Applied University (BAU) %
% Date of programming: 2020 %
% -------------------------------------------------
% This demo only implements a standard version of CapSA for a minimization problem 
% of a standard test function on MATLAB (R2018).
% -------------------------------------------------	
% Note:
% Due to the stochastic nature of meta-heuristc algorithms, 
% different runs may produce slightly different results.
%____________________________________________________________________________________

% This function initialize the first population of chameleons
function pos=init(searchAgents,dim,u,l)

for i=1:searchAgents 
    for j=1:dim
        pos(i,j)=round(l(j)-rand()*(l(j)-u(j)));
    end
end