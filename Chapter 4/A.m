x=linspace(0.99,1.01,101);
f=power(x,8)-8.*power(x,7)+28.*power(x,6)-56.*power(x,5)+70.*power(x,4)-56.*power(x,3)+28.*power(x,2)-8.*x+1;
g=(((((((x-8).*x+28).*x-56).*x+70).*x-56).*x+28).*x-8).*x+1;
h=power(x-1,8);
plot(x,f,'-r',x,g,'-b',x,h,'-k')
legend('f(x)','g(x)','h(x)')
title('Figure A')