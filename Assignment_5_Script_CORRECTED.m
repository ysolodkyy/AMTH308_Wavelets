clear all;
clc
close all;

%[h0,h1,h2,h3] = solve('h0^2+h1^2+h2^2+h3^2=1', 'h0+h1+h2+h3=sqrt(2)', 'h0-h1+h2-h3=0', 'h1-2*h2+3*h3=0')
%%
%{
syms x y

%ans1=x^2+y^2 +(1/sqrt(2)-x)^2 +(1/sqrt(2)-y)^2 ==1

%ans2 = 1/(2*sqrt(2))==x+y

figure()
ezplot(x^2+y^2 +(1/sqrt(2)-x)^2 +(1/sqrt(2)-y)^2 ==1, [-1/2,1])
hold on


ezplot(1/(2*sqrt(2))==x+y, [-1/2,1]) % should this be with positive slope? I don't trust this

hold off
% so to get the intersect values, I need to set the two eqs equal
clear x y
[x,y] = solve('x^2+y^2 +(1/sqrt(2)-x)^2 +(1/sqrt(2)-y)^2 =1','1/(2*sqrt(2))=x+y')

%}
%% plot the absolute value of the symbol 


[h0,h1,h2,h3] = solve('h0^2+h1^2+h2^2+h3^2=1', 'h0+h1+h2+h3=sqrt(2)', 'h0-h1+h2-h3=0', 'h1-2*h2+3*h3=0','h1-2*h2+3*h3=0')


syms z

h=[h0(2) h1(2) h2(2) h3(2)];


p=0;
for k = 0:3
p = p+1/sqrt(2)*h(k+1)*exp(-2*pi*j*z*k); % this may not be entirely correct. Consider doing symsubs
end

p

figure(1)
ezplot(abs(p), [0,1]) % this worked
title("Symbol P(z)");

%% PLot the scaling function see page 12

syms t

Phi_0 = heaviside(t)-heaviside(t-1);
figure(2); 
ezplot(Phi_0, [0,3])


for m=1:7
    Phi=0;
    for k=0:3
        Phi=Phi+sqrt(2)*h(k+1)*subs(Phi_0,t,2*t-k);
    end
figure(m+3)
ezplot(Phi, [0,3]);
title("Scaling Function Phi(t)");
Phi_0 = Phi;

%pause(.3);
end


%% mother wavelet see page 13

Psi = 0;

 for k=0:3
        Psi=Psi+(-1)^(k)*sqrt(2)*h(4-k)*subs(Phi,t,2*t-k);
 end

figure(m+3)
ezplot(Psi, [0,3]);
title("Mother Function Psi(t)");
