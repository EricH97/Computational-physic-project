function [x_nodes, w] = GaussLegendre(n)
%Computes the nodes and weights for Gauss-Legendre quadrature on the
%interval [-1, 1]
%Parameters
%       n = order of the quadrature
%Outputs
%       x_nodes = column vector containing the nodes (zeros of the Legendre
%       polynomial)
%       w = row vector containing the corresponding weights of the quadrature 

k=1:n;
beta=sqrt(k.^2 ./ (4 * k.^2 - 1));
alpha=zeros(size(beta)+[0,1]);

%Golub-Welsch algorithm
[Q,Lambda] = eig(diag(beta,1)+diag(alpha,0)+diag(beta,-1));
[x_nodes,i] = sort(diag(Lambda));
Qtop = Q(1,i);

w = 2*Qtop.^2;

end

