%define time scale
t0 = 0;
h = 0.001;
tf = 100;

%define parameters
a = [20, 1];

%define initial conditions
x0 = 0.2;
y0 = 0.2;
cond = [x0, y0];

%plot choice-- "0" for time courses /// "1" for trajectory 
plttgl = 0;

stl = '-';


rk4(t0, h, tf, cond, a, plttgl, 1, stl);


function rk4(t0, h, tf, y0, a, plttgl, stop, stl)

    function yout = orbit(t, y, a) 
    
    dy = zeros(2,1);  % column vector
    dy(1) = y(2);
    dy(2) = a(1)*(1-y(1)^2)*y(2)-y(1)+a(2)*sin(t);
    
    yout = dy;
   
    end

t=t0;

N=ceil(tf/h);

y1= zeros(N, 1);
y2= zeros(N, 1);


y1(1)=y0(1);
y2(1)=y0(2);

for i=1:N
    
    t(i+1)=t(i)+h;
    
    prntr = 100*(t(i+1)/tf);
    
    tm = t(i) + h/2;
    
    ys = [y1(i),              y2(i)];
    s1 = orbit(t(i),ys, a);
    
    ys2=[ys(1)+(h/2)*s1(1),   ys(2)+(h/2)*s1(2)];
    s2 = orbit(tm,ys2, a);
                                   
    ys3=[ys(1)+(h/2)*s2(1),   ys(2)+(h/2)*s2(2)];
    s3 = orbit(tm,ys3, a);
                                        
    ys4=[ys(1)+h*s3(1),       ys(2)+h*s3(2)];
    s4 = orbit(t(i+1),ys4, a);
    
    
    y1(i+1) = y1(i) + (h/6)*(s1(1) + 2*s2(1) + 2*s3(1) + s4(1));
    y2(i+1) = y2(i) + (h/6)*(s1(2) + 2*s2(2) + 2*s3(2) + s4(2));
    
end
    
figure(1);
if plttgl == 0 
   plot(t, y1);
   hold on
    if stop == 1
        hold off
        xlabel('t');
        ylabel('x(t)');
        title('Oscillations of x(t)--Functional Time Course');
    end

elseif plttgl == 1
    plot(y1, y2, stl)
    hold on
    if stop == 1
        hold off
        xlabel('x(t)');
        ylabel('y(t) = dx/dt');
        title('Phase Plane Trajectory--Stable Periodic Feedback Cyle');
    end
end
end

