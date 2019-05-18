function [F, a_coef, U] = force(Q, R, d, N, x_nodes, w)
%This function computes the force between two metal spheres with total
%charge of Q, a radius of R and a distance of d. We set 1/(4*pi*epsilon) =
%1.
%Parameters
%   Q = total charge on the spheres (same on both spheres)
%   R = radius of the spheres (same on both spheres)
%   d = distance between the speheres
%   N = The degree of the highest order Legendre polynomial taken into account 
%   when computing the surface charge density.
%   x_nodes = node points for the Gauss-Legendre quadrature
%   w = weights for the Gauss-Legendre quadrature.
%Outputs:
%   F = magnitude of the force between the spheres. 
%   a_coef = the coefficients of the Legendre polynomials related to the
%   surface charge density
%   U = potential energy between the spheres
%Note that we let epsilon = 1

A_matrix  = zeros(N, N);
b_vect = zeros(N, 1);
%The scaled integrals are stored to the matrix I_matrix
I_matrix = zeros(N+1, N+1); %Note that now the first row/column corresponds to 
%the Legendre polynomial n =0!!
I_matrix(1, 1) = 2*pi*R^3*2*(R/d);

%The first coefficient is dictated by the total charge
a0 = Q/(4*pi*R^2); 

for l = 1 : N
    for j = 1 : N
        %Let's evaluate the function to integrate at the node points
        fun_A_at_x_nodes = R^(j+1)/(2*j+1)*LP(l, x_nodes).*LP(j, -(d + R*x_nodes)./sqrt(d^2+R^2+2*d*R*x_nodes))./sqrt(d^2 + R^2 + 2*d*R*x_nodes).^(j+1);
        A_matrix(l, j) = w*fun_A_at_x_nodes;
        I_matrix(l+1,j+1) = 2*pi*R^3*A_matrix(l, j);
        if l == j
            b_vect(l) = -a0*2*(-1)^j/(2*j+1)*(R/d)^(j+1);
            A_matrix(l, l) = A_matrix(l, l) + 2/(2*j+1)^2; 
            I_matrix(l+1, 1) = 2*pi*R^3*2*(-1)^j/(2*j+1)*(R/d)^(j+1);
            I_matrix(1, l+1) = I_matrix(l+1, 1);
        end
    end
end


%Let's solve the resulting equation
a_coef = A_matrix\b_vect; %Column vector

%We prepend a0 as the first element
a_coef = [a0; a_coef];
a_coef(1);
I_matrix(1, 1);

%Computing the interaction potential energy between the spheres (doesn't
%include the self energies);
U = potEnergy(a_coef, I_matrix);


%Now, we can compute the force by fixing the charge distribution and
%monitor the interaction potential energy once the distance is slightly tweeked. 
%We apply the formula: f'(x) = 1/12h*(f(x-2*h) - 8*f(x-h) + 8*f(x+h) -f(x+2h))
h = 10^(-3);
delta_vect = [-2*h, -h, h, 2*h];
N_delta = length(delta_vect);
U_pot_vect = zeros(1, N_delta); %Potential at d-2h, d-h; d+h,d+2h

for k = 1 : N_delta
    dh = d + delta_vect(k);
    
    A_matrix  = zeros(N, N);
    b_vect = zeros(N, 1);
    %The scaled integrals are stored to the matrix I_matrix
    I_matrix = zeros(N+1, N+1); %Note that now the first row/column corresponds to 
    %the Legendre polynomial n =0!!
    I_matrix(1, 1) = 2*pi*R^3*2*(R/dh);
    for l = 1 : N
        for j = 1 : N
            %Let's evaluate the function to integrate at the node points
            fun_A_at_x_nodes = R^(j+1)/(2*j+1)*LP(l, x_nodes).*LP(j, -(dh + R*x_nodes)./sqrt(dh^2+R^2+2*dh*R*x_nodes))./sqrt(dh^2 + R^2 + 2*dh*R*x_nodes).^(j+1);
            A_matrix(l, j) = w*fun_A_at_x_nodes;
            I_matrix(l+1,j+1) = 2*pi*R^3*A_matrix(l, j);
            if l == j
                b_vect(l) = -a0*2*(-1)^j/(2*j+1)*(R/dh)^(j+1);
                A_matrix(l, l) = A_matrix(l, l) + 2/(2*j+1)^2; 
                I_matrix(l+1, 1) = 2*pi*R^3*2*(-1)^j/(2*j+1)*(R/dh)^(j+1);
                I_matrix(1, l+1) = I_matrix(l+1, 1);
            end
        end
    end
    %Compute the potential based on the tweeked integrals but keep the
    %surface charge density the same.
    U_pot_vect(k) = potEnergy(a_coef, I_matrix);
end
%Computing the force
F = (-1/(12*h)*(U_pot_vect(1) -8*U_pot_vect(2) + 8*U_pot_vect(3) - U_pot_vect(4)));

end

