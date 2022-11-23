x_1=[-3.5 , -3 , -2.5 , -2 , -1.75 , -1.5 , -1.25 , -1 , -0.875 , -0.75 , -0.625 , -0.5 , 0 , 0.5 , 0.625 , 0.75 , 0.875 , 1 , 1.25 , 1.5 , 1.75 , 2 , 2.5 , 3 , 3.5 ]';
y_1=zeros(25,1);
scatter(x_1,y_1,3,'filled','red');
hold on;
x_2=[-0.1875, -0.125, -0.0625, 0.0625, 0.125, 0.1875 ]';
y_2=zeros(6,1);
scatter(x_2,y_2,3,'filled','blue');
legend('Normal numbers','Subnormal numbers');
title('picture B');