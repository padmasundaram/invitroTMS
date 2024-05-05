function E = bemf4_surface_field_electric(c, Center, Area, eps0, SN, iintegralx, iintegraly, iintegralz)
%   Computes true continuous electric field on a surface facet due to
%   charges on ALL OTHER facets including accurate neighbor integrals.
%   Self-terms causing discontinuity may not be included
%
%   Copyright SNM 2017-2019    
   
    %  FMM 2019   
    %----------------------------------------------------------------------
    %   Fields plus potentials of surface charges (potential not used)
    %tic
    eps             = 1e-1;
    pg              = 2;
    srcinfo.sources = Center';
    srcinfo.charges = (c.*Area)';
    U               = lfmm3d(eps, srcinfo, pg);
    E               = -U.grad'/(4*pi);        
    %----------------------------------------------------------------------   
    %   Precise integration           
    %   Undo the effect of the m-th triangle charge on neighbors and
    %   add precise integration instead    
    %     for m = 1:M             
    %         index = ineighbor(m, :);          
    %         E(index, :)= E(index, :) + c(m)*[iintegralx(:, m) iintegraly(:, m) iintegralz(:, m)];
    %     end                         
    %   This is the vectorized implementation of the above code (Michael Tsuk, MathWorks, Inc.)
    cix = c.*iintegralx.';
    ciy = c.*iintegraly.';
    ciz = c.*iintegralz.';
    ci = [cix(:), ciy(:), ciz(:)];
    E = E + SN*ci;
    E = E/eps0;
    %toc 
end
