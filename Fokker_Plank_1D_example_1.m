%%% Sept 15, 2025 & IPT, MGL
%
function [] = Fokker_Plank_1D_example_1
%%%
%%% --- 1D Fokker-Plank equation with the following initial and boundary conditions
% du/dt = -d(A*u)/dx + d^2(B*u)/dx^2 with BC: u(0,t) = t
%                                             u(1,t) = 1 + t;
%                                         IC: u(x,0) = x, (0<x<1)
% A(x,u,t) = -1 & B(x,u,t) = 1.;
% with these A & B, one has du/dt = du/dx + d^2u/dx^2.
%
% analytic solution: u(x,t) = x + t, [1].
% The finite difference we use here is the backward-time centered space (BTCS). 
%%% ---
% matrix equation: 
% (I - lambda1*u_x - lambda2*u_xx)*u^(n+1) = u^(n) + 
%                                            lambda_1*b1^(n+1) + 
%                                            lambda_2*b2^(n+1)
%                                              
% A is tridiagonal matrix.
% (I - lambda1*u_x - lambda2*u_xx) is evolution matrix.
% b1^(n+1) = [-u^(n+1)_{a} 0,...0, u^(n+1)_{b}]^(T), (T = transpose). 
% b2^(n+1) = [u^(n+1)_{a} 0,...0, u^(n+1)_{b}]^(T), (T = transpose). 
%
% lambda1 = dt/(2*dx), lambda2 = dt/dx^2.
%
% Reference book: D. Bradie, A friendly introduction to numerical analysis 
% Reference [1]:   M. Dehghan and M. Tatari, Phys. Scr. 74 (2006) 310–316          
%
% Written by Tsogbayar Tsednee (PhD), Institute of Physics and Technology, Mongolian Academy of Sciences 
% Contact: tsog215@gmail.com
% Date: Sept 15, 2025
%
format short;
clc; clear; 
%
a = 0.; % x(0)
b = 1.0; % x(N+1)
N = 8;  % 
dx = (b-a)/N;
% grid and initial condition
% ---
x = zeros(N+1,1); % total number of points is N+1
for i = 1:N+1
    x(i) = a + (i-1)*dx;
end
%x;
% ---
ti = 0.; % t(0)
tf = 1.000; % t(N)
Nt = 32.; 
dt = (tf-ti)/Nt;
t = zeros(Nt,1); % total number of points is Nt
%
lambda1 = dt/(2.*dx);
lambda2 = dt/dx^2;
%
%  matrix equation is
% (I + lambda1*u_x - lambda2*u_xx)*u^(n+1) = u^(n) + 
%                                            lambda_1*b1^(n+1) + 
%                                            lambda_2*b2^(n+1)
%
u_mat_dx1 = zeros(N+1,N+1);
u_mat_dx2 = zeros(N+1,N+1);
for i = 2:N
    u_mat_dx2(i,i-1) = 1.0;
    u_mat_dx2(i,i) = -2.;
    u_mat_dx2(i,i+1) = 1.;
%
    u_mat_dx1(i,i-1) = -1.;
    u_mat_dx1(i,i) = 0.;
    u_mat_dx1(i,i+1) = 1.;
end
%u_mat_dx1;
%
u_mat_dx2(1,1) = -2.;
u_mat_dx2(N+1,N+1) = -2. ;
u_mat_dx2(1,2) = 1.;
u_mat_dx2(N+1,N) = 1.;
%
u_mat_dx2 = u_mat_dx2(1:N+1,1:N+1);
u_mat_dx2 = lambda2*u_mat_dx2 ;
%
u_mat_dx1(1,1) = 0.;
u_mat_dx1(N+1,N+1) = 0.;
u_mat_dx1(1,2) = 1.;
u_mat_dx1(N+1,N) = -1.;
%
u_mat_dx1 = u_mat_dx1(1:N+1,1:N+1) ;
u_mat_dx1 = lambda1*u_mat_dx1;
%
unit_I = eye(N+1);
% at t0 = 0
w0 = x; % % at t0 = 0.00 & initial
evolution_mat = unit_I - u_mat_dx1 - u_mat_dx2;
%%%
%b_n = zeros(N-1,1); 
w0 = w0(2:N); % taking into account the BC 
%
for j = 1:Nt
    t(j) = ti + (j-0)*dt;
    b_1n = [-t(j); zeros(N-3,1);  1. + t(j)];        
    b_2n = [ t(j); zeros(N-3,1);  1. + t(j)];    
%    
    w_old = w0 + lambda1.*b_1n + lambda2.*b_2n;
    w_new = evolution_mat(2:N,2:N)\(w_old  ); % % BTCS scheme
    w0 = w_new  ;
    %
    %figure(1)
    %hold on
    %plot(x(2:N), w0, 'r-'); drawnow;    
    %plot(x(2:N), x(2:N)+t(j), 'bo'); drawnow;
    %xlabel('$x$','interpreter','latex') % ,'fontsize',16
    %ylabel('$u(x,t)$','interpreter','latex') % ,'Rotation', 
    %hold off
    %
end
%%%
time = tf; % 
w_sol = [time; w_new; 1. + time];
u_exact = x + time;
%
figure(2)
hold on
plot(x, w_sol, 'b', LineWidth=1.5)
plot(x, u_exact, 'ro', LineWidth=1.5)
hold off
%axis([0. 1. .00 time])
set(gca,'FontSize',18)
xlabel('$x$','interpreter','latex') % ,'fontsize',16
ylabel('$u(x,t)$','interpreter','latex') % ,'Rotation', 
box on 
%%%
%rms(w_sol - u_exact);

return
end