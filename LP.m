function [LP_vect] = LP(n, x)
%(Standard) nth order Legendre polynomials on the interval x \in [-1, 1].
%Params
%       n = order of the legendre Polynomial (<= 10 !!!)
%       x = vector of x coordinates (we use x = cos theta)
%Output
%       LP_vect = nth order Legendre polynomial evaluated at the locations
%       x; Note that this is a column vector

[Nx, Ny] = size(x); %Number of x coordinates
x = reshape(x, [Nx, Ny]); %to column vector
LP_vect = zeros(Nx, Ny); %column vector

if n == 0
    LP_vect = ones(Nx, Ny);
elseif n == 1
    LP_vect = x;
elseif n == 2
    LP_vect = 0.5*(3*x.^2 - 1);
elseif n == 3
    LP_vect = 0.5*(5*x.^3 - 3*x);
elseif n == 4
    LP_vect = 1/8*(35*x.^4 - 30*x.^2 + 3);
elseif n == 5
    LP_vect = 1/8*(63*x.^5 - 70*x.^3 + 15*x);
elseif n == 6
    LP_vect = 1/16*(231*x.^6 - 315*x.^4 + 105*x.^2 -5);
elseif n == 7
    LP_vect = 1/16*(429*x.^7 - 693*x.^5 + 315*x.^3 - 35*x);
elseif n == 8
    LP_vect = 1/128*(6435*x.^8 - 12012*x.^6 + 6930*x.^4-1260*x.^2 + 35);
elseif n == 9
    LP_vect = 1/128*(12155*x.^9 - 25740*x.^7 + 18018*x.^5 - 4620 *x.^3 + 315*x);
elseif n == 10
    LP_vect =  1/256*(46189*x.^10 - 109395*x.^8 + 90090*x.^6 - 30030*x.^4 + 3465*x.^2 - 63);
end



end

