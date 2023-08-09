%% program do skrętu.
% Aksamit Michał
%
clear;
clc;
close all;

%% łuk 
a =- 100;

punkt1 = [-246.922; -35];
punkt2 = [-246.922; 35];

if a >= 0
    r1 = a + 201;
    r2 = r1 - 402;
    r3 = r1 - 75.5;
    r4 = r1 - 306.5;
    
    odchylenie = 45/(r1/200);
    odchylenie1 = 45/(r4/200);   
    
    rozdzielczosc = 1/(r1/50);
    
    if (r2<0) && (r2>-50)
        r4 =  sqrt( (r2)^2 + (95.5)^2);
        odchylenie1 = 45/(r4/135);
    elseif (r2 <= -50) && (r2>-150)
         r4 =  - sqrt( (r2)^2 + (95.5)^2) ;
        odchylenie1 = 45/(r4/200);
        punkt1 = [-246.922;  31];
        punkt2 = [-246.922; -31];
    elseif r2 <= -150 && r2 >= -201
         r4 = r2 + a/2;
        odchylenie1 = 45/(r4/200);
        punkt1 = [-246.922;  31];
        punkt2 = [-246.922; -31];
    end
    
    if r3 < 220
        r3 = r1 - (a/2);
    end
else
    r2 = a - 201;
    r1 = r2 + 402;
    r3 = r2 + 306.5;
    r4 = r2 + 75.5;
    
    odchylenie = 45/(r2/200);
    odchylenie1 = 45/(r3/200);
    
    rozdzielczosc = 1/(r2/250);
    
     if (r1>0) && (r1<50)
        r3 =  -sqrt( (r1)^2 + (95.5)^2);
        odchylenie1 = 45/(r3/135);
    elseif (r1 >= 50) && (r1 < 150)
         r3 = sqrt( (r1)^2 + (95.5)^2) ;
        odchylenie1 = 45/(r3/200);
        punkt1 = [-246.922;  31];
        punkt2 = [-246.922; -31];
    elseif r1 >= 150 && r1 < 201
         r3 = r1 - a/2;
        odchylenie1 = 45/(r3/200);
        punkt1 = [-246.922;  31];
        punkt2 = [-246.922; -31];
     end
    
     if r4 > -220
        r4 = r2 - (a/2);
     end   
end

prosta = sqrt((punkt2(1) - punkt1(1))^2 + (punkt2(2) - punkt1(2))^2);

if a>=0
    theta = 2*asind(prosta/(2*r1));
else
    theta = 2*asind(prosta/(2*r2));
end

angle = -(theta/2)+180: rozdzielczosc: (theta/2)+180+(theta/90);
x1 = r1*cosd(angle);
y1 = r1*sind(angle);

y1_real = -(x1 + a + 80);
x1_real = y1;

x2 = r2*cosd(angle);
y2 = r2*sind(angle);

y2_real = (x2 + a + 80);
x2_real = -y2;

if a >= 0
    angle1 = -(theta/2) + 180 + odchylenie: rozdzielczosc: (theta/2) + 180 + odchylenie + (theta/90);
    angle2 = -(theta/2) + 180 - odchylenie: rozdzielczosc: (theta/2) + 180 - odchylenie + (theta/90);
else
    angle1 = -(theta/2) + 180 - odchylenie1: rozdzielczosc: (theta/2) + 180 - odchylenie1 + (theta/90);
    angle2 = -(theta/2) + 180 + odchylenie1: rozdzielczosc: (theta/2) + 180 + odchylenie1 + (theta/90);
end

x3 = r3*cosd(angle1);
y3 = r3*sind(angle1);

y3_real = -(x3 + a + 80);
x3_real = y3;

x4 = r3*cosd(angle2);
y4 = r3*sind(angle2);

y4_real = -(x4 + a + 80);
x4_real = y4;

if a >= 0
    angle3 = -(theta/2) + odchylenie1 + 180: rozdzielczosc: (theta/2) + odchylenie1+(theta/90)+180;
    angle4 = -(theta/2) - odchylenie1 + 180: rozdzielczosc: (theta/2) - odchylenie1+(theta/90)+180;
else
    angle3 = -(theta/2) - odchylenie + 180: rozdzielczosc: (theta/2) - odchylenie +(theta/90)+180;
    angle4 = -(theta/2) + odchylenie + 180: rozdzielczosc: (theta/2) + odchylenie +(theta/90)+180; 
end

x5 = r4*cosd(angle3);
y5 = r4*sind(angle3);

y5_real = (x5 + a - 50);
x5_real = -y5;

x6 = r4*cosd(angle4);
y6 = r4*sind(angle4);

y6_real = (x6 + a - 50);
x6_real = -y6;

figure(1);
plot(x1, y1, '-o', 'color', 'red');
axis equal; % ustawienie takiej samej skali dla osi x i y
hold on;
plot(x2, y2, '-o', 'color', 'blue');
plot(x3, y3, '-o', 'color', 'green');
plot(x4, y4, '-o', 'color', 'green');
plot(x5, y5, '-o', 'color', 'yellow');
plot(x6, y6, '-o', 'color', 'yellow');
plot(0, 0, 'mo');
grid on;
hold off;

