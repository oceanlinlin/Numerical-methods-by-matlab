%------四阶经典龙格-库塔公式-----
function [x,y] = m4rkutta(df,xspan,y0,h)
%----四阶经典经典龙格-库塔公式解常微分方程 y'= f(x,y),y(x0)=y0----
%----格式：[x,y] = m4rkutta(df,xspan,y0,h),df为函数f(x,y)表达式，xspan为求解区间[x0,xn],y0为初值,h为步长，x为节点，y为数值解
format long
x = xspan(1):h:xspan(2);
y(1)=y0;
for n=1:(length(x)-1)
    k1 = feval(df,x(n),y(n));
    k2 = feval(df,x(n)+h/2,y(n)+h/2*k1);
    k3 = feval(df,x(n)+h/2,y(n)+h/2*k2);
    k4 = feval(df,x(n)+h,y(n)+h*k3);
    y(n+1) = y(n)+ h*(k1+2*k2+2*k3+k4)/6;
end

end
%调用格式如下
%df = @(x,y)3*y/(1+x);
%[x,y] = m4rkutta(df,[0 1],1,0.2)
