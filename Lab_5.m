clear all
clc
syms x1 x2;
fun_b=x1.^3+6*x1.^4-3*x2.^2+2*x2.^6
%x1.x2+3*x2+2(x1.^2)+x2.^2
GRADF=gradient(fun_b)
Hessf=hessian(fun_b)
[X1,X2]=meshgrid([-3:0.1:3],[-3:0.1:3]);
y=f(X1,X2);
contour(X1,X2,y)
figure
plot3(X1,X2,y)
%%
tol=0.001;
% Initial Guess (Choose Initial Guesses):
x1(1) = double(1);
x2(1) = double(1);
i=1;
d0 = -gradientf(x1(i),x2(i));
J=-d0;
dk=[d0];
beta=0;
B=eye(2);
Deltax=[5;5];
while ((norm(Deltax)>tol) && (i<200))

    dk(:,i)=-B*gradientf(x1(i),x2(i));

    S(1,i) =golden(0,1,0.001,x1(i),x2(i),dk(1,i));
    S(2,i) =golden(0,1,0.001,x1(i),x2(i),dk(2,i));
    x1(i+1) = double(x1(i)+S(1,i)*dk(1,i));
    x2(i+1) = double(x2(i)+S(2,i)*dk(2,i));

    Deltax=double([x1(i+1)-x1(i);x2(i+1)-x2(i)]);
    DeltaG=double(gradientf(x1(i+1),x2(i+1))-gradientf(x1(i),x2(i)));
    %david fletcher powell
    %B=double(B+(Deltax*Deltax')/(Deltax'*DeltaG)-(B*DeltaG*(B*DeltaG)')/(DeltaG'*B*DeltaG));
    %
    B=double(B+(DeltaG*DeltaG')/(DeltaG'*Deltax)-(B*Deltax*(B*Deltax)')/(Deltax'*B*Deltax));
    i = i+1;
end
raspn=double(vpa([x1;x2],4));
raspn(:,length(raspn)-1)
%%
figure
[X1,X2]=meshgrid([-3:0.1:3],[-3:0.1:3]);
contour(X1,X2,f(X1,X2))
hold on
for i=1:length(raspn)
    q=raspn(:,i);
    x1=vpa(q(1),4);
    x2=vpa(q(2),4);
    plot(q(1),q(2),'*')
    hold on
end
