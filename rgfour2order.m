%Gopesh Gaba 2020MCB1236
%RK4 method to solve IVP of 2nd order ODE with 1 independent and 1 dependent variable with
%y"=g(x,y,y')

clc
clear all

t0 = input('enter initial value of x\n');%initial value of x
u0 = input('enter initial value of u at initial x\n');%value of u at initial x
v0 = input('enter initial value of du/dx=v at initial x\n');%value of v=u' at initial x

syms u;
syms v;
syms x;%the above 3 lines initialize the variables

f(x,u,v)=v;%as v = u'
g(x,u,v)=input('enter dv/dx=d2u/dx2=g(x,u,v)\n');%to take input of v' or u"

h=input('enter step size\n');
n=input('enter number of steps\n');

F(x,u,v)=[f(x,u,v);g(x,u,v)];%column vector function to make coding easier

t=t0:h:t0+h*n;%row matrix containing points at which code will approximate the values of u and v
A=zeros(size(2,n+1));%matrix with 2 rows, first contains aproximations of u at corresponding t and second contains approximations of v at corresponding t
A(1,1)=u0;
A(2,1)=v0;

for i=1:n
    K1=vpa(h.*F(t(i),A(1,i),A(2,i)));
    K2=vpa(h.*F(t(i)+h/2,A(1,i)+K1(1,1)/2,A(2,i)+K1(2,1)/2));
    K3=vpa(h.*F(t(i)+h/2,A(1,i)+K2(1,1)/2,A(2,i)+K2(2,1)/2));
    K4=vpa(h.*F(t(i)+h,A(1,i)+K3(1,1),A(2,i)+K3(2,1)));
    A(:,i+1) = vpa(A(:,i)) + vpa((K1+2.*K2+2.*K3+K4)/6);
end%to compute approximations of u and u' at corresponding t

vpa(A);

tn = table(transpose(t),transpose(A(1,:)),transpose(A(2,:)));
tn = renamevars(tn,["Var1","Var2","Var3"],["x","u","du/dt or v"]);
tn%table to print the x and the corresponding u and u'

hold on
plot(t,A(1,:));
plot(t,A(2,:));
hold off% to plot u and v with respect to x
