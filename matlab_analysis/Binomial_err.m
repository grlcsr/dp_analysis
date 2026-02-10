%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [weight,p_err] = Binomial_err(x)
%
% simple function to compute the weigths and standard deviation
% associated to the input vector, presumed to have in each cell a
% Binomial distributed quantity.
% Weights are 1/sigma^2
%
% Possibly useful to properly weight the bin content in a histogram.
%
% WARNING: whenever the cell content is 0, the weight is ZERO and errors
% are also set to ZERO, to avoid fooling the graphical display!!!
%
% Author: Massimo Caccia, April 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [weight,p_err] = Binomial_err(x)
%   include later on a check on the dimension and apply it to 1-D only
%   vectors....
    n = length(x);
    Tot=sum(x);
    for j=1:n;
        if(x(j)==0);
            weight(j)=0;
            p_err(j)=0;
        else
            p=x(j)/Tot;
            p_err(j) = sqrt(Tot*p*(1-p));
            weight(j)=1/p_err(j)^2;
        end
        
    end
end

