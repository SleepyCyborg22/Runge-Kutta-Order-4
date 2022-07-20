%Gopesh Gaba
%2020MCB1236
%Runge Kutta Method order 4 for IVP of 1st order ODE with 1 dependent and 1 independent
%variable with y'=f(x,y)

clc
clear all

x0=input('enter initial x co-ordinate\n');%to have user input initial x coordinate
y0=input('enter initial y co-ordinate\n');%to have user input initial y coordinate

syms x;%to initialize variable x
syms y;%to initialize variable y

h=input('Enter step size\n');%to have user input the step size
n=input('Enter number of steps\n');%to have user input the number of steps

str=input('Enter f(x,y) suct that du/dt=f(t,u)\n');%to have user enter the function f(x,y) such that dy/dx=f(x,y) in form of a string

f = inline(str,'x','y');%to convert the string entered by user into a function

t = x0:h:x0+h*n;%to store the x coordniates.
u=zeros();%to store the y coordinates
u(1)=y0;%intial y coordinate



for i=1:n%Loop that runs n times to calculate the next n y coordinates
    K1=h*f(t(i),u(i));%to calculate K1
    K2=h*f(t(i)+h/2,u(i)+K1/2);%to calculate K2
    K3=h*f(t(i)+h/2,u(i)+K2/2);%to calculate K3
    K4=h*f(t(i)+h,u(i)+K3);%to calculate K4
    u(i+1)=u(i)+(K1+2*K2+2*K3+K4)/6;%to find u(i+1) using u(i),K1,K2,K3,K4
end%to compute the values of K1, K2, K3, K4(i.e.) the weighted average of slopes and use them to calculate the next y coordinate.


tn = table(transpose(t),transpose(u));
tn = renamevars(tn,["Var1","Var2"],["X Coordinate","Y Coordinate"]);
tn%table to print the x and the corresponding y coordniate

plot(t,u)%to plot t on x axis and u on y axis
