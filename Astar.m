%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Astar Application
clear;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define the map and initial values
% 0 for initial path
% 1 for open path
% 2 for close path
% 10 for final path

% define map

Map = [2, 2, 2, 2, 2, 2, 2;
       2, 1, 1, 1, 1, 1, 2;
       2, 1, 2, 2, 1, 2, 2;
       2, 1, 2, 0, 2, 2, 2;
       2, 1, 1, 2, 2, 10, 2;
       2, 1, 1, 1, 1, 1, 2;
       2, 2, 2, 2, 2, 2, 2];
   

% define initial and final coordinate
[xi,yi] = find(Map == 0);
[xf,yf] = find(Map == 10);
initial = [xi,yi];
final = [xf,yf];

var_coor = initial;
Map(initial(1,1),initial(1,2)) = 2;

%% define the map and initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% algortihm

step_matrix = [-1,-1 ; % Left Up
                -1,0 ; % Up
                -1,1 ; % Right Up
                0,-1 ; % Left               
                 0,1 ; % Right 
                1,-1 ; % Left Down                
                 1,0 ; % Down 
                 1,1]; % Right Down 
             
path_matrix = initial;
temp_path_matrix = initial;
j=1;
t=0;
while (var_coor(1,1) ~= final(1,1) || var_coor(1,2) ~= final(1,2))
  
    s=0;
    
    distance_matrix = [1000;
                       1000;
                       1000;
                       1000;
                       1000;
                       1000;
                       1000;
                       1000]; 
                   
    for i = 1 : 1 : length(step_matrix(:,1))
        
        int_step = var_coor + step_matrix(i,:);
        
        s = s + Map(int_step(1,1),int_step(1,2));
        
        if(Map(int_step(1,1),int_step(1,2)) == 10)
            distance_matrix(i,1) = (int_step(1,1)-final(1,1))^2+(int_step(1,2)-final(1,2))^2;
            break;
        end
        
        if(s == 16)
            Map(var_coor(1,1),var_coor(1,2)) = 2;
                        
            if(var_coor(1,1) == initial(1,1) && var_coor(1,2) == initial(1,2))
                j=j+1; 
                var_coor = temp_path_matrix(j,:);
                                
                t=t+1;
                if(t > 1000)
                   error('No Way!');
                end
                
            else
            var_coor = initial;            
            
            end
            path_matrix = initial;
        end
    
        
        if (Map(int_step(1,1),int_step(1,2)) ~= 2 )
            distance_matrix(i,1) = (int_step(1,1)-final(1,1))^2+(int_step(1,2)-final(1,2))^2;
            
        end
    end
    
        
    minimum = min(distance_matrix); 
    
    
    for i = 1 : 1 : length(distance_matrix(:,1))
        
        if(minimum == 1000)
            break;
            
        elseif (minimum ==  distance_matrix(i,1))            
            var_coor = var_coor + step_matrix(i,:);            
            Map(var_coor(1,1),var_coor(1,2)) = 2;
            break;
            
        end
        
    end
    
    path_matrix = [path_matrix;var_coor];
    temp_path_matrix = [temp_path_matrix;var_coor];
        
end
path_matrix
%% algortihm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Astar Application
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
