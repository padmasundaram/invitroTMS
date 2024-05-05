function [integralxc, integralyc, integralzc, integralpc,...
          integralxe, integralye, integralze, integralpe,...
          EC, PC] = meshneighborints(P, t, normals, Area, Center, RnumberE, ineighborE, numThreads)
%   Accurate integration for electric field/electric potential on neighbor facets
%   Copyright SNM 2017-2020
    tic 
    N = size(t, 1);
    integralxe      = zeros(RnumberE, N);    %   exact Ex integrals for array of neighbor triangles 
    integralye      = zeros(RnumberE, N);    %   exact Ey integrals for array of neighbor triangles 
    integralze      = zeros(RnumberE, N);    %   exact Ez integrals for array of neighbor triangles 
    integralpe      = zeros(RnumberE, N);    %   exact potential integrals for array of neighbor triangles 

    integralxc      = zeros(RnumberE, N);    %   center-point Ex integrals for array of neighbor triangles 
    integralyc      = zeros(RnumberE, N);    %   center-point Ey integrals for array of neighbor triangles 
    integralzc      = zeros(RnumberE, N);    %   center-point Ez integrals for array of neighbor triangles 
    integralpc      = zeros(RnumberE, N);    %   center-point potential integrals for array of neighbor triangles 

    gauss       = 25;   %   number of integration points in the Gaussian quadrature  
                        %   for the outer potential integrals
                        %   Numbers 1, 4, 7, 13, 25 are permitted 
    %   Gaussian weights for analytical integration (for the outer integral)
    if gauss == 1;  [coeffS, weightsS, IndexS]  = tri(1, 1); end;
    if gauss == 4;  [coeffS, weightsS, IndexS]  = tri(4, 3); end;
    if gauss == 7;  [coeffS, weightsS, IndexS]  = tri(7, 5); end;
    if gauss == 13; [coeffS, weightsS, IndexS]  = tri(13, 7); end;
    if gauss == 25; [coeffS, weightsS, IndexS]  = tri(25, 10); end;
    W           = repmat(weightsS', 1, 3);

    %   Main loop for analytical double integrals (parallel, 24 workers)
    %   This is the loop over columns of the system matrix
    tic
    %parpool(16);
    parpool(numThreads);
    parfor n = 1:N                  %   inner integral; (n =1 - first column of the system matrix, etc.)        
        r1      = P(t(n, 1), :);    %   [1x3]
        r2      = P(t(n, 2), :);    %   [1x3]
        r3      = P(t(n, 3), :);    %   [1x3]  
        index   = ineighborE(:, n); %   those are non-zero rows of the system matrix for given n
        ObsPoints   = zeros(RnumberE*IndexS, 3);    %   to compute RnumberE outer integrals numerically
        I           = zeros(RnumberE, 3);           %   for rhe field
        IP          = zeros(RnumberE, 1);           %   for the potential
        %   Accurate electric-field integrals
        for q = 1:RnumberE
            num = index(q);
            for p = 1:IndexS
                ObsPoints(p+(q-1)*IndexS, :)  = coeffS(1, p)*P(t(num, 1), :) +  coeffS(2, p)*P(t(num, 2), :) +  coeffS(3, p)*P(t(num, 3), :);
            end
        end    
        J = potint2(r1, r2, r3, normals(n, :), ObsPoints); %   Outer integral is computed analytically, for all inner IntPoints  
        for q = 1:RnumberE
            I(q, :) = sum(W.*J([1:IndexS]+(q-1)*IndexS, :), 1); 
        end
        I(1, :)          = 0;           %   self integrals will give zero
        integralxe(:, n) = -I(:, 1);    %   accurate integrals, entries of non-zero rows of n-th column
        integralye(:, n) = -I(:, 2);    %   accurate integrals, entries of non-zero rows of n-th column
        integralze(:, n) = -I(:, 3);    %   accurate integrals, entries of non-zero rows of n-th column  
        %   Center-point electric-field integrals
        temp    = repmat(Center(n, :), RnumberE, 1) - Center(index, :); %   these are distances to the observation/target triangle
        DIST    = sqrt(dot(temp, temp, 2));                             %   single column                
        I       = Area(n)*temp./repmat(DIST.^3, 1, 3);                  %   center-point integral, standard format    
        I(1, :) = 0;                                                    %   self integrals will give zero
        integralxc(:, n) = -I(:, 1);    %   center-point integrals, entries of non-zero rows of n-th column
        integralyc(:, n) = -I(:, 2);    %   center-point integrals, entries of non-zero rows of n-th column
        integralzc(:, n) = -I(:, 3);    %   center-point integrals, entries of non-zero rows of n-th column        
        %   Accurate electric-potential integrals
        [JP, dummy] = potint(r1, r2, r3, normals(n, :), ObsPoints);     %   JP was calculated without the area Area(n)      
        for q = 1:RnumberE      
            IP(q)   = sum(W(:, 1).*JP([1:IndexS]+(q-1)*IndexS), 1);
        end
        integralpe(:, n) = +IP;                     %   accurate integrals (here is without the area!)        
        %   Center-point electric-potential integrals         
        IPC     = Area(n)./DIST;                                        %   center-point integral, standard format 
        IPC(1)  = 0;                                                    %   this must be zero
        integralpc(:, n) = +IPC;                                        %   center-point approximation
    end 
    delete(gcp('nocreate'));
    
    %%  Define useful sparse matrices EC, PC (for GMRES speed up)    
    N               = size(t, 1);
    const           = 1/(4*pi);  
    integralc       = zeros(RnumberE, N);    %   normal integral component for array of neighbor triangles (center point) - to speed up GMRES
    integrale       = zeros(RnumberE, N);    %   normal integral component for array of neighbor triangles (exact) - to speed up GMRES
    for n = 1:N                  %   inner integral; (n =1 - first column of the system matrix, etc.)             
        index = ineighborE(:, n); %   those are non-zero rows of the system matrix for given n
        integrale(:, n)  =       +(integralxe(:, n).*normals(index, 1) + ...
                                   integralye(:, n).*normals(index, 2) + ...
                                   integralze(:, n).*normals(index, 3)); 
        integralc(:, n)  =       +(integralxc(:, n).*normals(index, 1) + ...
                                   integralyc(:, n).*normals(index, 2) + ...
                                   integralzc(:, n).*normals(index, 3));
    end 
    ii  = ineighborE;
    jj  = repmat([1:N], RnumberE, 1);
    EC  = sparse(ii, jj, const*(-integralc + integrale));               %   almost symmetric
    PC  = sparse(ii, jj, const*(-integralpc + integralpe));             %   almost symmetric   
    
    MainLoopParallelTime = toc
end