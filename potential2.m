function [V_matrix] = potential2(x1, x2, R, a_coef, x_matrix, y_matrix)
%Potential around AND INSIDE the spheres at the locationcs xvect, y_vect. We set
%1/(4*pi*epsilon) = 1;
%Parameters
%   x1 = location of the first sphere
%   x2 = location of the 2nd sphere
%   R = radius of the spheres
%   a_coef = coefficients of the Legendre polynomials
%   x_matrix = x coordinates of the places at which we want to compute the
%   potential
%   y_matrix = y coordinates of the places at which we want to compute the
%   potential
%Outputs
%   V_matrix = potential at points x_vect, y_vect

[nrows, ncols] = size(x_matrix);
V_matrix = zeros(nrows, ncols);

%Highest order of Legendre polynomials
N_L = length(a_coef) - 1;

%Distances to the centers of the spheres
r1_matrix = sqrt((x_matrix - x1).^2 + y_matrix.^2);
r2_matrix = sqrt((x_matrix - x2).^2 + y_matrix.^2);

idx_inside1 = find(r1_matrix < R);
idx_inside2 = find(r2_matrix < R); %We assume that spheres don't overlap.
idx_outside = find(r1_matrix >= R & r2_matrix >= R);

%Cosines of the angles 
cos1_matrix = -(x_matrix - x1)./r1_matrix;
cos2_matrix = (x_matrix - x2)./r2_matrix;


for i = 0 : N_L
    %4*pi is due to (4*pi*eps)^(-1) = 1
    %Potential for points outside the spheres
    V_matrix(idx_outside) = V_matrix(idx_outside) + 4*pi*a_coef(i+1)*R^(i + 2)/(2*i + 1)*(LP(i, cos1_matrix(idx_outside))./r1_matrix(idx_outside).^(i + 1) + LP(i, cos2_matrix(idx_outside))./r2_matrix(idx_outside).^(i + 1));
    %Potential for points inside the sphere 1
    V_matrix(idx_inside1) = V_matrix(idx_inside1) + 4*pi*a_coef(i+1)/(2*i + 1)*(R^(-i + 1)*LP(i, cos1_matrix(idx_inside1)).*r1_matrix(idx_inside1).^i  + R^(i + 2)*LP(i, cos2_matrix(idx_inside1))./r2_matrix(idx_inside1).^(i + 1));
    %Potential for points inside the sphere 2
    V_matrix(idx_inside2) = V_matrix(idx_inside2) + 4*pi*a_coef(i+1)/(2*i + 1)*(R^(i + 2)*LP(i, cos1_matrix(idx_inside2))./r1_matrix(idx_inside2).^(i+1)  + R^(-i + 1)*LP(i, cos2_matrix(idx_inside2)).*r2_matrix(idx_inside2).^(i));
end

end

