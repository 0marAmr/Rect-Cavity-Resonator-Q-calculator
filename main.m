function main
clc;
clear;
close all;
format compact;

%cavity resonator dimensions in cm
global a; 
global b; 
global d;

%other cavity parameteers
global epsilon_r;            %relative permitivity 
global sigma;                %metal conductivity
global Q_die;                %loss tangent at specific freq.

% epsilon_r = 2.5;
% sigma = 5.8*10^7;
% loss_tanget = 0.00004;
% Q_die = 1/loss_tanget;
% a = 0.05;
% b = 0.04;
% d = 0.07;

%get inputs from user
a = input("a (in cm)= ")*10^-2;
b = input("b (in cm)= ")*10^-2;
d = input("d (in cm)= ")*10^-2;
epsilon_r = input("relative permeability ε_r = ");
sigma = input("metal conductivity σ = ");
Q_die = 1/input("loss_tanget tan(δ) = ");

%constants
global c_0;
global epsilon_0;
global u_0;

c_0 = 3*10^8;
epsilon_0 = 8.8542* 10^-12;
u_0 = 4*pi*10^-7;

f_res_dominant =  resonance_freq(1,0,1); %#ok<*NOPRT,*NASGU>
Q_dominant = Q_total_10l(1,2*pi*f_res_dominant)

r = sqrt( (3*(a^2)*(b^2))/(a^2-b^2));

if d > r % next mode is TE102
   f_res_next = resonance_freq(1,0,2);
   Q_next = Q_total_10l(2,2*pi*f_res_next);
else
   f_res_next = resonance_freq(0,1,1);
   Q_next = Q_total_01l(1,2*pi*f_res_next);
end

print_output(f_res_dominant,Q_dominant,f_res_next,Q_next,r);
end

function [f_resonance] = resonance_freq(m,n,l)

global a; 
global b; 
global d;
global epsilon_r; 

c_0= 3*10^8;
c = c_0/sqrt(epsilon_r);
kx = m*pi/(a);
ky = n*pi/(b);
beta = l*pi/(d);
f_resonance = (c/(2*pi)) * sqrt(kx^2 + ky^2 +beta^2); 
end

%TE101 , TE102
function [Q_10l] = Q_total_10l(l,w_res)
global a; 
global b; 
global d;
global epsilon_r; 
global epsilon_0;
global u_0;
global sigma;  
global Q_die;

k_c= pi/a;
R_s = sqrt( (w_res * u_0)/(2*sigma));
Q_num = w_res *(((epsilon_r*epsilon_0*a*b*d)/2) * ((w_res*u_0)/k_c)^2);
Q_den = R_s * (2*(((a*l)/d)^2)*a*b + 2*d*b + a*d + (((a*l)/d)^2)*a*d);

Q_cond_10l = Q_num/Q_den;
dummy = 1/Q_cond_10l + 1/Q_die;
Q_10l = 1/dummy;
end

%TE011
function [Q_01l] = Q_total_01l(l,w_res)

global a; 
global b; 
global d;
global epsilon_r; 
global epsilon_0;
global u_0;
global sigma;
global Q_die;

R_s=sqrt( (w_res * u_0)/(2*sigma));
k_c=pi/b;

Q_num = w_res *(((epsilon_r*epsilon_0*a*b*d)/2) * ((w_res*u_0)/k_c)^2);
Q_den = R_s * (2*(((b*l)/d)^2)*a*b + 2*a*d + b*d + (((b*l)/d)^2)*b*d);
Q_cond_01l = Q_num/Q_den;

dummy = 1/Q_cond_01l + 1/Q_die;
Q_01l = 1/dummy;
end

function print_output(f_res_dominant,Q_dominant,f_res_next,Q_next,r)
global d;
if d >r
    next_mode = "TE102";
else
    next_mode = "TE011";
end

BW_dominant = 1/Q_dominant;
BW_next = 1/Q_next;

fprintf("\n==================================\n");
fprintf("The dominant mode: TE101\n");
fprintf("the 1st resonance frequency = %f GHz\n",f_res_dominant*10^-9);
fprintf("the quality factor Q (TE101) = %f \n",Q_dominant);
fprintf("the fractional bandwidth (TE101) = %e \n",BW_dominant);
fprintf("==================================\n");
fprintf("The next mode: %s\n",next_mode);
fprintf("the 2nd resonance frequency = %f GHz\n",f_res_next*10^-9);
fprintf("the quality factor Q(%s) = %f \n",next_mode,Q_next);
fprintf("the fractional bandwidth (%s) = %e \n",next_mode,BW_next);
fprintf("==================================\n");
end