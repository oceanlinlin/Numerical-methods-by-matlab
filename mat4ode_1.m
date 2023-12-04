%-----https://ww2.mathworks.cn/help/matlab/math/solve-bvp-with-unknown-parameter.html

%用bvp4c,即有限差分方法求解出二阶常微分方程的边值问题
%y''-4x^2y=1+2x,x∈(0,2),y'(0)=1,y(2)=2


%使用区间为[0,2]的10点网格，初始估计值函数调用bvpinit
solinit = bvpinit(linspace(0,2,10),@guess);

%使用ODE函数、边界条件函数和初始估计值调用bvp4c
sol = bvp4c(@mat4ode,@mat4bc,solinit);

%对解进行绘图
% 使用 deval 计算 bvp4c 在区间 [0,2] 中的 11 个点处计算的解
xint = linspace(0,2,11);
Sxint = deval(sol,xint);
%对两个解分量进行绘图。
plot(xint,Sxint)
axis([0 2 -4 4])
title('Eigenfunction of Mathieu''s Equation.') 
xlabel('x')
ylabel('y')
legend('y','y''')

%dydx = mat4ode(x,y,lambda);%x是自变量,y是因变量,lambda为特征值的未知参数
%----用代换法y1=y,y2=y'将二阶常微分方程组方程转化为一阶方程组
function dxdy = mat4ode(x,y)
dxdy = [y(2) 
        4*x^2*y(1)+(1+2*x)];
end

%编写边界条件代码
function res = mat4bc(ya,yb)
res = [ya(2)-1
    %yb(2)
    yb(1)-2
    ];
end

%创建初始估计值,可以估计对应的y值为指数函数，也可以估计为线性函数x。验证过了，这两种结果相同
function yinit = guess(x)
yinit = [exp(x)
         exp(x)];
end
