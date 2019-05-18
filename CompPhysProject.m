clear all
close all

%Note that initially the code was quite slow. By defining an own
%function for the Legendre polynomials, using Gauss-Legendre quadrature 
%instead of the integral function and exploiting the parallelism in built
%in the vector operations, I managed to get the code run about 200 faster!!

%Note that througout the code 1/(4*pi*epsilon_0) = 1

%Let's use the first 11 Legendre polynomials in approximating the charge
%distribution
R = 1; %Radius of the spheres
Q = 1; %Total charge of the spheres.
a0 = Q/(4*pi*R^2); 
N = 10; %in addition to P_0, we use N other Legendre polynomials
A_matrix = zeros(N, N);
b_vect = zeros(N, 1);

%The distance between the spheres
d =2;


%We precompute the nodes and weights of the Gauss-Legendre quadrature
[x_nodes, w] = GaussLegendre(100); %order of 100 should be enough
%We compare the results for different N
N_vect = [2,3,5,7,10]; 
N_theta = 100;
sigma_matrix = zeros(length(N_vect), N_theta); %Matrix for storing the surface charge densities
tic
for idx_N = 1 : length(N_vect)
    N = N_vect(idx_N);
    A_matrix = zeros(N, N);
    b_vect = zeros(N, 1);
    tic
    for l = 1 : N
        for j = 1 : N
            %Let's evaluate the function to integrate at the node points
            fun_A_at_x_nodes = R^(j+1)/(2*j+1)*LP(l, x_nodes).*LP(j, -(d + R*x_nodes)./sqrt(d^2+R^2+2*d*R*x_nodes))./sqrt(d^2 + R^2 + 2*d*R*x_nodes).^(j+1);
            A_matrix(l, j) = w*fun_A_at_x_nodes;
            if l == j
                b_vect(l) = -a0*2*(-1)^j/(2*j+1)*(R/d)^(j+1);
                A_matrix(l, l) = A_matrix(l, l) + 2/(2*j+1)^2; 
            end
        end
    end
    disp('Numerical integration')
    toc
    
    %Let's solve the resulting matrix equation
    tic 
    a_coef = A_matrix\b_vect; %Column vector
    disp('Solving coefs')
    toc
    %Note that the numerical integrations takes about 100 times longer than
    %solving the matrix equation.
    
    
    %We prepend a0 as the first element
    a_coef = [a0; a_coef];
    %The corresponding surface charge density as a function of theta
    N_theta = 100;
    theta_vect = linspace(0, pi, N_theta);
    LP_matrix = zeros(N_theta, N+1);
    for i = 0 : N
        LP_matrix(:, i+1) = LP(i, cos(theta_vect))'; %to column vector
    end
    sigma_matrix(idx_N, :) = (LP_matrix *a_coef)';
end


%Plotting the surface charge density for different N when d = 3.
figure;
hold on
for idx = 1 : length(N_vect)
    plot(theta_vect, 4*pi*sigma_matrix(idx, :), 'Linewidth', 2.5)
end
set(gca, 'ticklength', 2*get(gca, 'ticklength'), 'Linewidth', 1.3, 'Fontsize', 14)
xlabel('$\theta$ (rad)', 'Interpreter', 'latex', 'Fontsize', 22)
ylabel('$4\pi \sigma(\theta)$', 'Interpreter', 'latex', 'Fontsize', 22)
%title('Surface charge density','Interpreter', 'latex', 'Fontsize', 20)
leg = legend(['$N =$ ', num2str(N_vect(1))], ['$N =$ ', num2str(N_vect(2))], ['$N =$ ', num2str(N_vect(3))],['$N =$ ', num2str(N_vect(4))], ['$N =$ ', num2str(N_vect(5))],  'Location', 'southwest');
set(leg, 'Interpreter', 'latex', 'Fontsize', 20)
box on 
grid on
xlim([0, pi])
%ylim([0, 1.5])
hold off

%Let's make a nice 3D plot of the surface charge density. N = 10
N_sph = 50;
[Y,Z,X] = sphere(N_sph) ; %// generate coordinates of a sphere    

[azimuth,elevation,r] = cart2sph(Y,Z,X) ;  %// Convert cartesian coordinates to spherical referential
LP_matrix = zeros(N_sph+1, N+1);
for i = 0 : N
    LP_matrix(:, i+1) = LP(i, cos(pi/2 - elevation(:, 1))); %column vector
end
sigma3D = LP_matrix *a_coef;
sigma3D = repmat(sigma3D, 1, N_sph+1);
 %// surface charge density
figure;
hold on
surf(X + d/2, Y, Z, 4*pi*sigma3D)
surf(X - d/2, Y, Z, flipud(4*pi*sigma3D));
axis equal           ; %// set the axis ratio so the sphere appear as a sphere
shading interp       ; %// small refinement to not see the grid
set(gca, 'ticklength', 2*get(gca, 'ticklength'), 'Linewidth', 1.3, 'Fontsize', 14)
xlabel('$x$', 'Interpreter','LaTex','FontSize', 20);
ylabel('$y$', 'Interpreter','LaTex','FontSize', 20);
zlabel('$z$', 'Interpreter','LaTex','FontSize', 20);
h = colorbar;
ylabel(h, 'Surface charge density $4\pi\sigma$',  'Interpreter','LaTex','FontSize', 20) 
colormap(jet)
box on 
%Let's check whether the total potential is constant on the surface of the sphere
%1.

x1 = -d/2;
x2 = d/2;
x_vect = (x1 + R*cos(theta_vect))';
y_vect =  R*sin(theta_vect)';
[V_vect] = potential2(x1, x2, R, a_coef, x_vect, y_vect);

%Yes, it is (roughly) constant!! -> Now, all the equations are correct.
figure;
hold on
plot(theta_vect, V_vect, 'r', 'Linewidth', 2)
set(gca, 'ticklength', 2*get(gca, 'ticklength'), 'Linewidth', 1.3, 'Fontsize', 14)
xlabel('$\theta_1$ (rad)', 'Interpreter', 'latex', 'Fontsize', 22)
ylabel('Potential $V$ on $\partial S_1$', 'Interpreter', 'latex', 'Fontsize', 22)
%title('Potential on sphere 1','Interpreter', 'latex', 'Fontsize', 20)
%leg = legend(['$n_\mathrm{highest} =$ ', num2str(N2-1)], ['$n_\mathrm{highest} =$ ', num2str(N)],  'Location', 'northeast');
%set(leg, 'Interpreter', 'latex', 'Fontsize', 20)
box on 
grid on
xlim([0, pi])
hold off


%Let's compute the potential around the spheres
N_pot = 500;
x_vect = linspace(-3, 3, N_pot);
y_vect = linspace(-3, 3, N_pot);
[X, Y] = meshgrid(x_vect, y_vect);
[V_matrix] = potential2(x1, x2, R, a_coef, X, Y);

figure;
hold on
surf(X, Y, V_matrix)
axis tight;
shading flat;
view([0 0 90]);
set(gca,'FontSize',14, 'Linewidth', 1.3);
xlabel('$x$', 'Interpreter','LaTex','FontSize', 20);
ylabel('$y$', 'Interpreter','LaTex','FontSize', 20);
h = colorbar;
ylabel(h, 'Potential $V$',  'Interpreter','LaTex','FontSize', 20) 
%title('Potential around the spheres', 'Interpreter','LaTex','FontSize', 20) 
box on
grid on
colormap(jet)



%Let's determine the force between the charged balls as a function of
%distance
%Let's also determine the surface charge density at the closest and
%furthest point as a function of distance
%For reference, we compute the force based on the image charge method
N_d = 100;
d_vect = linspace(2.01*R, 15*R, N_d);
F_vect = zeros(1, N_d);
F_image = zeros(1, N_d);
U_pot_vect = zeros(1, N_d);
surf_charge_0 = zeros(1, N_d); %Furthest point theta = 0
surf_charge_pi = zeros(1, N_d); %Closest point theta = pi
%We take into account Legendre polynomials up to order 10.
N = 10;
LP_vect_0 = zeros(1, N+1);
LP_vect_pi = zeros(1, N+1);
for i = 0 : N
    LP_vect_0(:, i+1) = LP(i, cos(0)); 
    LP_vect_pi(:, i+1) = LP(i, cos(pi)); 
end


for k = 1 : N_d
    d = d_vect(k);
    [F, a_coef, U] = force(Q, R, d, N, x_nodes, w);
    F_vect(k) = F;
    surf_charge_0(k) = LP_vect_0*a_coef;
    surf_charge_pi(k) = LP_vect_pi*a_coef;
    U_pot_vect(k) = U;
    F_image(k) = forceFromImages(Q, R, d, 1000);
end

%Plotting the surface charge density as a function of distance
figure;
hold on
plot(d_vect, 4*pi*surf_charge_0, 'r', 'Linewidth', 2.5)
plot(d_vect, 4*pi*surf_charge_pi, 'b', 'Linewidth', 2.5)
set(gca, 'ticklength', 2*get(gca, 'ticklength'), 'Linewidth', 1.3, 'Fontsize', 14)
xlabel('Distance $d$', 'Interpreter', 'latex', 'Fontsize', 22)
ylabel('$4\pi \sigma$', 'Interpreter', 'latex', 'Fontsize', 22)
%title('Force as a function of distance','Interpreter', 'latex', 'Fontsize', 20)
leg = legend('$\theta = 0$ rad', '$\theta = \pi$ rad',  'Location', 'northeast');
set(leg, 'Interpreter', 'latex', 'Fontsize', 17)
box on 
grid on
xlim([2, 15])
ylim([0, 1.5])
hold off


%Plotting the force as a function of distance
figure;
hold on
plot(d_vect, F_vect, 'r', 'Linewidth', 2.5)
plot(d_vect, F_image, 'k--', 'Linewidth', 2.5)
plot(d_vect, 1./(d_vect.^2), 'b--', 'Linewidth', 2.5)
plot([2], 1./(2.^2)*0.6149, 'og', 'Linewidth', 2.5)
set(gca, 'ticklength', 2*get(gca, 'ticklength'), 'Linewidth', 1.3, 'Fontsize', 14)
xlabel('Distance $d$', 'Interpreter', 'latex', 'Fontsize', 22)
ylabel('Force $F$', 'Interpreter', 'latex', 'Fontsize', 22)
%title('Force as a function of distance','Interpreter', 'latex', 'Fontsize', 20)
leg = legend('Variational approach', 'Method of images', 'Point charge approx.', 'Result by Thomson',  'Location', 'northeast');
set(leg, 'Interpreter', 'latex', 'Fontsize', 17)
box on 
grid on
hold off


%Potential energy as a function of distance
figure;
hold on
plot(d_vect, U_pot_vect, 'r', 'Linewidth', 2.5)
plot(d_vect, 1./(d_vect), 'b--', 'Linewidth', 2.5)
set(gca, 'ticklength', 2*get(gca, 'ticklength'), 'Linewidth', 1.3, 'Fontsize', 14)
xlabel('$d$', 'Interpreter', 'latex', 'Fontsize', 22)
ylabel('$\epsilon U /Q^2$', 'Interpreter', 'latex', 'Fontsize', 22)
title('Potential as a function of distance','Interpreter', 'latex', 'Fontsize', 20)
leg = legend('Potential', '$\frac{1}{d}$',  'Location', 'northeast');
set(leg, 'Interpreter', 'latex', 'Fontsize', 20)
box on 
grid on
hold off

%LogLog scale
%Plotting the force as a function of distance
figure;
loglog(d_vect, F_vect, 'r', 'Linewidth', 2.5)
hold on
loglog(d_vect, F_image, 'k--', 'Linewidth', 2.5)
loglog(d_vect, 1./(d_vect.^2), 'b--', 'Linewidth', 2.5)
loglog([2], 1./(2.^2)*0.6149, 'og', 'Linewidth', 2.5)
set(gca, 'ticklength', 2*get(gca, 'ticklength'), 'Linewidth', 1.3, 'Fontsize', 14)
xlabel('Distance $d$', 'Interpreter', 'latex', 'Fontsize', 22)
ylabel('Force $F$', 'Interpreter', 'latex', 'Fontsize', 22)
%title('Force as a function of distance','Interpreter', 'latex', 'Fontsize', 20)
leg = legend('Variational approach', 'Method of images', 'Point charge approx.', 'Result by Thomson',  'Location', 'northeast');
set(leg, 'Interpreter', 'latex', 'Fontsize', 17)
box on 
grid on
xlim([2, 15])
hold off

%The corresponding surface charge density as a function of theta
tic
N_theta = 1000;
theta_vect = linspace(0, pi, N_theta);
LP_matrix = zeros(N_theta, N+1);
for i = 0 : N
    LP_matrix(:, i+1) = LP(i, cos(theta_vect))'; %to column vector
end
sigma = LP_matrix *a_coef;
sigma2 = LP_matrix(:, 1:3)*a_coef(1:3);%Surface charge density with only legendre polynomials n<=2.
toc

%Plotting the surface charge density
figure;
hold on
plot(theta_vect, sigma2, 'r', 'Linewidth', 2.5)
plot(theta_vect, sigma, 'b--', 'Linewidth', 2.5)
set(gca, 'ticklength', 2*get(gca, 'ticklength'), 'Linewidth', 1.3, 'Fontsize', 14)
xlabel('$\theta$', 'Interpreter', 'latex', 'Fontsize', 22)
ylabel('$\sigma(\theta)$', 'Interpreter', 'latex', 'Fontsize', 22)
title('Surface charge density, $d$=15','Interpreter', 'latex', 'Fontsize', 20)
leg = legend('3 polynomials', [num2str(N +1) ' polynomials'],  'Location', 'northeast');
set(leg, 'Interpreter', 'latex', 'Fontsize', 20)
box on 
grid on
xlim([0, pi])
hold off


%Let's determine the relative error in the force for different choices of N
%as a function of the distance between the centers of the spheres.
%We compare the result to the method of images that should be quite
%accurate if we use many enough image charges.

N_d = 30;
%Different N's to try
N_vect = [2, 3, 5, 7, 10];
d_vect = linspace(2.1*R, 15*R, N_d);
F_mat = zeros(length(N_vect), N_d);
rel_error_mat = zeros(length(N_vect), N_d); %Relative error of the force
F_image = zeros(1, N_d);
for idx_N = 1 : length(N_vect)
    N = N_vect(idx_N); %Highest order of the Legendre polynomials taken into account.
    for k = 1 : N_d
        d = d_vect(k); %Distance between spheres.
        %Force based on the variational approach
        [F, a_coef, U] = force(Q, R, d, N, x_nodes, w);
        F_mat(idx_N, k) = F;
        %Force based on the method of image charges.
        %This we use as a reference
        F_img = forceFromImages(Q, R, d, 1000);
        F_image(k) = F_img; %We do this N_d times. Doen't matter
        %Computing the relative error
        rel_error_mat(idx_N, k) = abs(F - F_img)/F_img;
    end
end

%Plotting the relative error for different value of N as a function of
%distance
%Plotting the surface charge density for different N when d = 3.
figure;
hold on
for idx = 1 : length(N_vect)
    plot(d_vect, rel_error_mat(idx, :), 'Linewidth', 2.5)
end
set(gca, 'ticklength', 2*get(gca, 'ticklength'), 'Linewidth', 1.3, 'Fontsize', 14)
xlabel('Distance $d$', 'Interpreter', 'latex', 'Fontsize', 22)
ylabel('Relative error $\Delta F/F_\mathrm{image}$', 'Interpreter', 'latex', 'Fontsize', 22)
%title('Surface charge density','Interpreter', 'latex', 'Fontsize', 20)
leg = legend(['$N =$ ', num2str(N_vect(1))], ['$N =$ ', num2str(N_vect(2))], ['$N =$ ', num2str(N_vect(3))],['$N =$ ', num2str(N_vect(4))], ['$N =$ ', num2str(N_vect(5))],  'Location', 'northeast');
set(leg, 'Interpreter', 'latex', 'Fontsize', 16)
set(gca, 'YScale', 'log') %semilogy
box on 
grid on
xlim([2, 15])
hold off

