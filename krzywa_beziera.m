%% program do krzywych beziera.
% Aksamit Michał
%
clear;
clc;
close all;

%% wstep 

% wsp. t : 0,1
syms t;

% liczba punktow
n = 7;

n_silnia = 1;
i_silnia = 1;
roznica_silnia = 1;
dwumian = 0;
punkt_Beziera  = 0;

krzywa_Beziera_x = 0;
krzywa_Beziera_y = 0;
krzywa_Beziera_z = 0;

Start = [-192; 75; -200];
Koniec = [-292; 75; -200];

if ( Start(1,1) - Koniec(1,1) ) == 0
    a = 0;
else
    a = atand( (Start(2,1) - Koniec(2,1)) / ( Start(1,1) - Koniec(1,1) ) );    
end    
    
if ( Start(1,1) < Koniec(1,1))
    P_1 = [ Start(1,1) - (20* cosd(a)); Start(2,1) - (26* sind(a)); Start(3,1)];
    P_2 = [ Start(1,1) - (26* cosd(a)); Start(2,1) - (26* sind(a)); Start(3,1) + 20];
    P_3 = [(Start(1,1) + Koniec(1,1))/2; (Start(2,1) + Koniec(2,1))/2; ((Start(3,1) + Koniec(3,1))/2)+100];
    P_4 = [ Koniec(1,1) + (20* cosd(a)); Koniec(2,1) + (26* sind(a)); Koniec(3,1) + 20];
    P_5 = [ Koniec(1,1) + (26* cosd(a)); Koniec(2,1) + (26* sind(a)); Koniec(3,1)];
else
    P_1 = [ Start(1,1) + (20* cosd(a)); Start(2,1) + (26* sind(a)); Start(3,1)];
    P_2 = [ Start(1,1) + (26* cosd(a)); Start(2,1) + (26* sind(a)); Start(3,1) + 20];
    P_3 = [(Start(1,1) + Koniec(1,1))/2; (Start(2,1) + Koniec(2,1))/2; ((Start(3,1) + Koniec(3,1))/2)+100];
    P_4 = [ Koniec(1,1) - (20* cosd(a)); Koniec(2,1) - (26* sind(a)); Koniec(3,1) + 20];
    P_5 = [ Koniec(1,1) - (26* cosd(a)); Koniec(2,1) - (26* sind(a)); Koniec(3,1)];
end

Punkt_x = [Start(1,1) P_1(1,1) P_2(1,1) P_3(1,1) P_4(1,1) P_5(1,1) Koniec(1,1)];
Punkt_y = [Start(2,1) P_1(2,1) P_2(2,1) P_3(2,1) P_4(2,1) P_5(2,1) Koniec(2,1)];
Punkt_z = [Start(3,1) P_1(3,1) P_2(3,1) P_3(3,1) P_4(3,1) P_5(3,1) Koniec(3,1)];

%% dwumian newtona

%{
for i=1 : n-1
    n_silnia = n_silnia * i;
end


for i = 0 : n-1
    for j = 0 : i 
        i_silnia = i_silnia * j;
        if i_silnia == 0
            i_silnia = 1;
        end   
    end
    
    for k = 0 : n-1-i 
        roznica_silnia = roznica_silnia * k;
        if roznica_silnia == 0
            roznica_silnia = 1;
        end   
    end
    
    dwumian = n_silnia/ (i_silnia * roznica_silnia);
    
   krzywa_Beziera_x = krzywa_Beziera_x + (dwumian * ((1-t)^(n-1-i))*(t^i)*Punkt_x(i+1));
   krzywa_Beziera_y = krzywa_Beziera_y + (dwumian * ((1-t)^(n-1-i))*(t^i)*Punkt_y(i+1));
   krzywa_Beziera_z = krzywa_Beziera_z + (dwumian * ((1-t)^(n-1-i))*(t^i)*Punkt_z(i+1));
end   

disp(krzywa_Beziera_x);
disp(krzywa_Beziera_y);
disp(krzywa_Beziera_z);
%}


for i=1 : n-1
    n_silnia = n_silnia * i;
end


figure (1);
hold on;
for s = 0:0.05:1 
    if ( Start(1,1) - Koniec(1,1) ) == 0
        break;
    end
    for i = 0 : n-1
        for j = 0 : i 
            i_silnia = i_silnia * j;
            if i_silnia == 0
                i_silnia = 1;
            end   
        end

        for k = 0 : n-1-i 
            roznica_silnia = roznica_silnia * k;
            
            if roznica_silnia == 0
                roznica_silnia = 1;
            end   
        end

       dwumian = n_silnia/ (i_silnia * roznica_silnia);

       krzywa_Beziera_x = krzywa_Beziera_x + (dwumian * ((1-s)^(n-1-i))*(s^i)*Punkt_x(i+1));
       krzywa_Beziera_y = krzywa_Beziera_y + (dwumian * ((1-s)^(n-1-i))*(s^i)*Punkt_y(i+1));
       krzywa_Beziera_z = krzywa_Beziera_z + (dwumian * ((1-s)^(n-1-i))*(s^i)*Punkt_z(i+1));
       
    end   
        
    plot3(krzywa_Beziera_x , krzywa_Beziera_y  ,krzywa_Beziera_z , '-o', 'color', 'blue');
    
    krzywa_Beziera_x = 0;
    krzywa_Beziera_y = 0;
    krzywa_Beziera_z = 0;
end
grid on;
hold off;

%{
figure (1);
hold on;
for s = 0:0.05:1 
  punkt_Beziera_x = (6479314710149097*s*(s - 1)^5)/17592186044416 - (1456084128273059*(s - 1)^6)/35184372088832 - (2370877802582691*s^5*(6*s - 6))/35184372088832 - (35563167038740365*s^2*(s - 1)^4)/35184372088832 + (32396573550745485*s^4*(s - 1)^2)/35184372088832 + (1456084128273059*s^6)/35184372088832;
  punkt_Beziera_y = 131*(s - 1)^6 - 786*s*(s - 1)^5 - 131*s^5*(6*s - 6) + 1965*s^2*(s - 1)^4 - 2620*s^3*(s - 1)^3 + 1965*s^4*(s - 1)^2 + 131*s^6;
  punkt_Beziera_z = 960*s*(s - 1)^5 - 160*(s - 1)^6 + 160*s^5*(6*s - 6) - 2100*s^2*(s - 1)^4 + 1200*s^3*(s - 1)^3 - 2100*s^4*(s - 1)^2 - 160*s^6;
  
  plot3(punkt_Beziera_x, punkt_Beziera_y ,punkt_Beziera_z,'-o','color', 'blue');
end
grid on;
hold off;
%}