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
% different runs may produce slightly different results.____________________________

function [fFitness,CapuchinPositions,cg_curve]=CapuchinAlg(noP,maxite,LB,UB,dim,fobj)
warning off; format long;

% % % CapSA main program
% f1 =  figure (1);
% set(gcf,'color','w');
% hold on
% xlabel('Iteration','interpreter','latex','FontName','Times','fontsize',10)
% ylabel('fitness value','interpreter','latex','FontName','Times','fontsize',10); 
% grid;
% %  
cg_curve=zeros(1,maxite);

%% % CapSA initialization

%Initialize the positions of Capuchins in the space
CapuchinPositions=init(noP,dim,UB,LB);

v=0.1*CapuchinPositions;% initial velocity
v0=zeros(noP,dim); 

% Calculate the fitness of initialCapuchins

    for i=1:noP
     CapuchinFitness(i,1)=fobj(CapuchinPositions(i,:));
    end
 
fitness = CapuchinFitness;
% Initial fitness of the random positions

[fFitness,index]=min(CapuchinFitness);

CapuchinBestPosition = CapuchinPositions; % Best position initialization 
gFoodPosition = CapuchinPositions(index,:); % initial global position 
%% % CapSA Parameters
bf=0.90;%Balance factor
cr=19.0;  %Modulus of elasticity
g=9.81;

% CapSA velocity updates
a1=1.250; a2=1.5;   

beta=[0.05 11 2];
wmax=0.9;
wmin=0.1;
%% % CapSA Main loop
    
for t = 1 : maxite

    % Life time convergence

        tau = beta(1) * exp(-beta(2) * ((t)/maxite)^beta(3));
        w   = wmax-(wmax-wmin)*(t/maxite)^1; 
 
        % CapSA velocity update
     
    for i=1:noP 
        for j=1:dim
        v(i,j)=  w* v(i,j) + ...
             1*  a1*(CapuchinBestPosition(i,j)- CapuchinPositions(i,j))*rand + ...
              1* a2*(gFoodPosition(j) - CapuchinPositions(i,j))*rand; 
 
        end        
    end
     
% CapSA position update

for i=1:noP
   if i<noP/2
      for j=1:dim
          if (rand()>=0.1)
               r=rand;
              if r<=0.15
                 CapuchinPositions(i,j) =  gFoodPosition(j)+    1*bf*((v(i,j).^2)*sin(2*rand()*1.5))/g;  % Jumping (Projection)
              elseif   r>0.15 && r<=0.30  
                  CapuchinPositions(i,j) =  gFoodPosition(j)+    1*cr*bf*((v(i,j).^2)*sin(2*rand()*1.5))/g;  % Jumping (Land)  
              elseif   r>0.30 && r<=0.9      
                  CapuchinPositions(i,j) =    CapuchinPositions(i,j) +  v(i,j); % Movement on the ground    
              elseif  r>0.9 && r<=0.95 
                 CapuchinPositions(i,j) =      gFoodPosition(j)  +  bf*1.0*sin(1*rand()*1.5);   % Swing % Local search  
              elseif   r>0.95 
                CapuchinPositions(i,j) =       gFoodPosition(j)  +  bf*1.0*(v(i,j)- v0(i,j));    % Climbing   % Local search
              end
        else
            CapuchinPositions(i,j) =             tau*(LB(j)  + rand *(UB(j)    - LB(j))); 
        end      
    end
% Let the followers follow the leaders (update their positions)
elseif i>=noP/2 && i<=noP 
        alphaPos=CapuchinPositions(i-1,:);
        followerPos=CapuchinPositions(i,:);
        CapuchinPositions(i,:)=(alphaPos+followerPos)/2; 
   end
end 
v0 = v;

for i=1:noP % relocation (Update, exploration)
        u=CapuchinPositions(i,:)>UB;
        l=CapuchinPositions(i,:)<LB;
     
         CapuchinPositions(i,:)= (CapuchinPositions(i,:).*~xor(u,l))+UB.*u+LB.*l;
    
         CapuchinFitness(i,1)=fobj (CapuchinPositions(i,:)) ;
         
            if CapuchinFitness(i,1)<fitness(i)
                CapuchinBestPosition(i,:)=CapuchinPositions(i,:);
                fitness(i)=CapuchinFitness(i,1);
            end 
end
%% Evaluate the new positions

[fmin,index]=min(fitness); % finding out the best positions  

% Updating gPosition and best fitness
if fmin < fFitness
    gFoodPosition = CapuchinBestPosition(index,:); % Update the global best positions
    fFitness = fmin;
end

% % %     % Display the iterative results
% 
%      outmsg = ['Iteration# ', num2str(t) , '  Fitness= ' , num2str(fFitness)];
%         disp(outmsg);
% % %    
% 
%    % Obtain the convergence curve
     cg_curve(t)=fFitness;
%      
%  if t>2
%      set(0, 'CurrentFigure', f1)
% 
%         line([t-1 t], [cg_curve(t-1) cg_curve(t)],'Color','b'); 
%         xlabel('Iteration');
%         ylabel('Best score obtained so far');
%         drawnow 
%  end
   
end

fFitness =fobj(gFoodPosition);
end