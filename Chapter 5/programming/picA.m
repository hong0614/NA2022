x = linspace(0,10,21);
y = [2.9 2.7 4.8 5.3 7.1 7.6 7.7 7.6 9.4 9.0 9.6 10.0 ...
    10.2 9.7 8.3 8.4 9.0 8.3 6.6 6.7 4.1];
plot(x,y,'o')
hold on;
grid on;
[a,s] = dls(x,y,3);

X = linspace(0,10,100001);
f = a(1) + a(2)*X + a(3)*X.*X;
plot(X,f,LineWidth=1.5);
hold on;

title("DLSA");
legend('discrete data','best approximation')