x = linspace(0,30);
y1 = 6.67+1.77167.*(x-0)+0.457833.*(x-0).*(x-6)-0.124778.*(x-0).*(x-6).*(x-10)+0.013566.*(x-0).*(x-6).*(x-10).*(x-13)-0.000978085.*(x-0).*(x-6).*(x-10).*(x-13).*(x-17)+4.1477e-05.*(x-0).*(x-6).*(x-10).*(x-13).*(x-17).*(x-20);
y2 = 6.67+1.57167.*(x-0)-0.0871667.*(x-0).*(x-6)-0.0152729.*(x-0).*(x-6).*(x-10)+0.00257908.*(x-0).*(x-6).*(x-10).*(x-13)-0.000204804.*(x-0).*(x-6).*(x-10).*(x-13).*(x-17)+8.6768e-06.*(x-0).*(x-6).*(x-10).*(x-13).*(x-17).*(x-20);
plot(x,y1,x,y2)