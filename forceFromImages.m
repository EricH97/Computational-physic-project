function [F] = forceFromImages(Q, R, d, N_i)
%We compute the force between two charged metallic spheres by applying the
%method of images. We set 1/(4*pi*epsilon) = 1
%   Q = total charge of the spheres
%   R = radius of the spheres
%   d = distance between the centers of the spheres
%   N_i = number of image charges taken into account
%   (>=2)
%Outputs
%   F = force between the spheres (epsilon = 1)

%Let's store the values of the image charges and their locations (with
%respect to the middle point between the sphere)
q_vect = zeros(1, N_i);
x_vect = zeros(1, N_i);

%The first image charge is located at the center of the spheres and we set
%its charge to be 1. Later the sum of image charges is scaled to equal the
%total charge
q_vect(1) = 1;
x_vect(1) = d/2;

for i = 2 : N_i
    q_vect(i) = -R/(d/2 + x_vect(i-1))*q_vect(i-1);
    x_vect(i) = d/2 - R^2/(d/2 + x_vect(i-1));
end

%Now, we scale the sum of charges so that it matches the total charge
q_vect = q_vect*Q/sum(q_vect);

%Now, we can compute the force between the image charges
F = 0;
for i = 1 : N_i
    for j = 1 : N_i
        F = F + q_vect(i)*q_vect(j)/(x_vect(i) + x_vect(j))^2;
    end
end
F = F;

end

