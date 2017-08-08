%INPUT 1----define temporal parameters
t0 = 0;
dt = 0.001;
tf = 50.0;
T = [t0, dt, tf];

%INPUT 2----define parameters
a= 10.0;
b= 24.4;
c= 8.0/3.0;
A = [a, b, c];

%INPUT 3----define initial conditions
x0 = 0.0;
y0 = 0.1;
z0 = 10.0;
init_cond = [x0, y0, z0];

%INPUT 4----choose plot output- 0 for oscillations// 1 for trajectories
plot = 1;


rk4(T, init_cond, A, plot, 1);


function rk4 = rk4(T, y0, A, plt, stop)
t0 = T(1);
h = T(2);
tf = T(3);

    function yout = orbit(t, y, A)
       
    a=A(1);
    b=A(2);
    c=A(3);
     
    dy = zeros(3,1);
    
    %INPUT 5---Define the x derivative wrt t // Using parameters a, b, c, d
    %as well as y(1) as x and y(2) as y
    dy(1) = a*(y(2)-y(1));
    
    % Define the y deritative wrt t similarly
    dy(2) = (b+sin(0.01*t))*y(1) - y(2) - y(1)*y(3);
    
    dy(3) = y(1)*y(2) - c*y(3);
    
    yout = dy;
   
    end

t=t0;

N=ceil(tf/h);

y1= zeros(N, 1);
y2= zeros(N, 1);
y3= zeros(N, 1);


y1(1)=y0(1);
y2(1)=y0(2);
y3(1)=y0(3);




for i=1:N
    
    t(i+1)=t(i)+h;
    
    prntr = abs(t(i+1)-t0)/abs(tf-t(i+1));
    
    tm = t(i) + h/2;
    
    
    ys = [y1(i),              y2(i), y3(i)];
                                        s1 = orbit(t(i),ys, A);
    
    
    
    ys2=[ys(1)+(h/2)*s1(1),   ys(2)+(h/2)*s1(2), ys(3)+(h/2)*s1(3)];
                                        s2 = orbit(tm,ys2, A);
                                        
                                        
    
    ys3=[ys(1)+(h/2)*s2(1),   ys(2)+(h/2)*s2(2), ys(3)+(h/2)*s2(3)];
                                        s3 = orbit(tm,ys3, A);
                                        
                                        
    
    ys4=[ys(1)+h*s3(1),       ys(2)+h*s3(2), ys(3)+h*s3(3)];
                                        s4 = orbit(t(i+1),ys4, A);
    
    
    y1(i+1) = y1(i) + (h/6)*(s1(1) + 2*s2(1) + 2*s3(1) + s4(1));
    y2(i+1) = y2(i) + (h/6)*(s1(2) + 2*s2(2) + 2*s3(2) + s4(2));
    y3(i+1) = y3(i) + (h/6)*(s1(3) + 2*s2(3) + 2*s3(3) + s4(3));
    
end

figure(1);
if plt == 0 
   plot(t, y1);
   hold on
    if stop == 1
        hold off
        xlabel('t');
        ylabel('x(t)');
        title('Oscillations of x(t)--Functional Time Course');
    end

elseif plt == 1
    plot(y1, y3, '-')
    hold on
    if stop == 1
        hold off
        xlabel('x');
        ylabel('z');
        title('Lorenz Trajectory: x(t) vs. z(t)');
    end
end
end