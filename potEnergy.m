function [U] = potEnergy(a_coef, I_matrix)
%This function computes the potential energy of the other charged sphere in
%the electric field created by the other sphere. This potential energy 
%describing the interactions can be used to compute the force between the
%spheres. We set 1/(4*pi*epsilon) = 1.
%   a_coef = column vector containing the weights of the Legendre
%   Polynomials related to the surface charge density.
%   I_matrix = matrix containing the integrals needed to compute the A
%   matrix. Note that these need to be scaled with 2*pi*R^3 BEFORE passing
%   to this function

U = 4*pi*a_coef'*I_matrix*a_coef;

end

