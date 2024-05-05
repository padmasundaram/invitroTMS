function [output] = mygrid(limit1, limit2, N)
    output = limit1 + [0:N-1]*(limit2-limit1)/N + 0.5*(limit2-limit1)/N;
end

