function [K,E] = ellipfast(m,kind)
%[K,E] = ellipfast(m,kind)
%   Computes the complete elliptic integral of the first kind K(m) and 
%   second kind E(m) for 0 <= m <= 1 with an error of less than 2e-8.
%   'kind' is either 1 or 2 or 'K' or 'E'.

% Source: Abramowitz 17.3.34 (page 591)

% S. H. Muller, 20/01/2003

if nargin > 1
    switch kind
    case {'k','K',1}
        K = ellipfastmex(m,1);
        return;
    case {'e','E',2}
        K = ellipfastmex(m,2);
        return;
    end
end

K = ellipfastmex(m,1);
E = ellipfastmex(m,2);