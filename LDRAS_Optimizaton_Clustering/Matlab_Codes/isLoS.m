%% FUNCTION TO KNOW IF THE RECEIVER IS IN LoS
function LoS = isLoS(j,k,mt,nt,x_p,y_p,x_bs,y_bs,wall_p_x,wall_p_y,slope,ordi,normal_to_x)

    LoS=1;

    % Find crossing walls
    if(j < mt)
        % East
        wall_sighted_in_x = ((wall_p_x < x_bs & wall_p_x > x_p));        
    elseif(j > mt)
        % West    
        wall_sighted_in_x = ((wall_p_x > x_bs & wall_p_x < x_p));        
    else
        % Same x-coordinate
        wall_sighted_in_x = (wall_p_x == x_bs);
    end
    
    if(k < nt)
        % North
        wall_sighted_in_y = ((wall_p_y < y_bs & wall_p_y > y_p));        
    elseif(k > nt)
        % South
        wall_sighted_in_y = ((wall_p_y > y_bs & wall_p_y < y_p));        
    else
        % Same y-coordinate
        wall_sighted_in_y = (wall_p_y == y_bs);
    end    
    
    for w=1:length(wall_p_x)
        if wall_sighted_in_x(w)
            x_wall_x = wall_p_x(w);
            y_wall_x = slope*x_wall_x + ordi;            
            build_y_1 = (wall_p_y(1:2:end) <= y_wall_x);     
            build_y_2 = (wall_p_y(2:2:end) >= y_wall_x);
            build_y = (build_y_1 & build_y_2);
            if sum(build_y)~=0
                LoS = 0;
                break;
            end                
        end
    end
    for w=1:length(wall_p_y)
        if wall_sighted_in_y(w)
            y_wall_y = wall_p_y(w);
            if ~normal_to_x
                x_wall_y = (y_wall_y - ordi)/slope;
            else
                x_wall_y = x_p;
            end            
            build_x_1 = (wall_p_x(1:2:end) <= x_wall_y);     
            build_x_2 = (wall_p_x(2:2:end) >= x_wall_y);
            build_x = (build_x_1 & build_x_2);            
            if sum(build_x)~=0
                LoS = 0;
                break;
            end            
        end
    end
end