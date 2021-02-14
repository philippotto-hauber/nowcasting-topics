function ystar = f_interpol(y)
% function to linearly interpolate missings in an NxT matrix 
ystar = y; 
for i = 1:size(y, 1)
    t = 2;
    while t < size(y, 2)
      ind_nextobs = 1;
        if isnan(ystar(i, t))  
            t_NaN = true;
            while t_NaN
                if isnan(ystar(i, t+ind_nextobs))
                    ind_nextobs = ind_nextobs + 1;
                else
                    t_NaN = false;
                end
            end
            ystar(i, t:t+ind_nextobs - 1) = ystar(i, t-1) + (1:ind_nextobs) * ((ystar(i, t + ind_nextobs) - ystar(i, t-1))/ind_nextobs);
        end
        t = t + ind_nextobs;
    end
end

% replace missings in period 1 and T with values from 2 or T-1
if any(isnan(y(:, 1)))
    ind_nan = isnan(y(:, 1));
    ystar(ind_nan, 1) = ystar(ind_nan, 2);
end

if any(isnan(y(:, end)))
    ind_nan = isnan(y(:, end));
    ystar(ind_nan, end) = ystar(ind_nan, end-1);
end

end