%% program do kinematyki hexapoda na prace inzynierska - wersja numeryczna.
% Aksamit Michał
%
clear;
clc;
close all;

%% dane wejsciowe
poz_konc = [-191;85.5;-207];

%starting position by angles
theta_1 = 135 - 90; 
theta_2 = 72 - 90;
theta_3 = -72;

% dimensions lengths of robot leg
L1 = 51;
L2 = 70;
L3 = 207;

% angle in iteration during time
actual_angle = [theta_1 theta_2 theta_3];

% pt and pg
difer_x =1;
difer_y =1;
difer_z =1;

% count iteration
k = 1;

% set weight for least squares
weight = [1; 2; 3];

%regularisation
R = [0;0;0];

% set bound for servos reasons
angleupperbound = [80;60;-20];
anglelowerbound = [-80;-60;-180];

% macierz punktów

   punkty = zeros(10000,3);
   punkty_B = zeros(10000,3);
   punkty_C = zeros(10000,3);

   B = zeros(3,3,10000);

%
  

%% petla glowna programu

 while ( sqrt(difer_x^2) > 0.001 || difer_x == 0 ) || ( sqrt(difer_y^2) > 0.001 || difer_y == 0) || ( sqrt(difer_z^2) > 0.001 || difer_z == 0)
            
        x_pos = L3*(sind(actual_angle(1,1))*sind(actual_angle(1,2))*sind(actual_angle(1,3)) - cosd(actual_angle(1,2))*cosd(actual_angle(1,3))*sind(actual_angle(1,1))) - L1*sind(actual_angle(1,1)) - L2*cosd(actual_angle(1,2))*sind(actual_angle(1,1));
        y_pos = L1*cosd(actual_angle(1,1)) - L3*(cosd(actual_angle(1,1))*sind(actual_angle(1,2))*sind(actual_angle(1,3)) - cosd(actual_angle(1,1))*cosd(actual_angle(1,2))*cosd(actual_angle(1,3))) + L2*cosd(actual_angle(1,1))*cosd(actual_angle(1,2));
        z_pos = L3*(cosd(actual_angle(1,2))*sind(actual_angle(1,3)) + cosd(actual_angle(1,3))*sind(actual_angle(1,2))) + L2*sind(actual_angle(1,2));
         
        %kinematyka prosta punktu B (przegub)
        xB_pos  = - L1*sind(actual_angle(1,1)) - L2*cosd(actual_angle(1,2))*sind(actual_angle(1,1));
        yB_pos = L1*cosd(actual_angle(1,1)) + L2*cosd(actual_angle(1,2))*cosd(actual_angle(1,2));
        zB_pos = L2*sind(actual_angle(1,2));

        %kinematyka prosta punktu C (przegub)
        xC_pos  = -L1*sind(actual_angle(1,1));
        yC_pos = L1*cosd(actual_angle(1,1));
        zC_pos = 0;

        efect_pos = [x_pos;y_pos;z_pos];

        %zebranie punktow do odpowiadajacych im tabel
        punkty(k,1) = x_pos;
        punkty(k,2) = y_pos;
        punkty(k,3) = z_pos;
        
        punkty_B(k,1) = xB_pos;
        punkty_B(k,2) = yB_pos;
        punkty_B(k,3) = zB_pos;
        
        punkty_C(k,1) =  xC_pos;
        punkty_C(k,2) =  yC_pos;
        punkty_C(k,3) =  zC_pos;

        % budowa regularyzacji
        R(1,1) = 1.0 * actual_angle(1,1)^2;
        R(1,2) = 1.0 * actual_angle(1,2)^2;
        R(1,3) = 1.0 * actual_angle(1,3)^2;

        %tworzenie jakobiana do obliczeń numerycznych
        J_num = zeros(3);

        J_num(1,1) = L3*(cosd(actual_angle(1,1))*sind(actual_angle(1,2))*sind(actual_angle(1,3)) - cosd(actual_angle(1,1))*cosd(actual_angle(1,2))*cosd(actual_angle(1,3))) - L1*cosd(actual_angle(1,1)) - L2*cosd(actual_angle(1,1))*cosd(actual_angle(1,2));
        J_num(2,1) = L3*(sind(actual_angle(1,1))*sind(actual_angle(1,2))*sind(actual_angle(1,3)) - cosd(actual_angle(1,2))*cosd(actual_angle(1,3))*sind(actual_angle(1,1))) - L1*sind(actual_angle(1,1)) - L2*cosd(actual_angle(1,2))*sind(actual_angle(1,1));
        J_num(3,1) = 1;
        J_num(1,2) = L3*(cosd(actual_angle(1,2))*sind(actual_angle(1,1))*sind(actual_angle(1,3)) + cosd(actual_angle(1,3))*sind(actual_angle(1,1))*sind(actual_angle(1,2))) + L2*sind(actual_angle(1,1))*sind(actual_angle(1,2));
        J_num(2,2) = - L3*(cosd(actual_angle(1,1))*cosd(actual_angle(1,2))*sind(actual_angle(1,3)) + cosd(actual_angle(1,1))*cosd(actual_angle(1,3))*sind(actual_angle(1,2))) - L2*cosd(actual_angle(1,1))*sind(actual_angle(1,2));
        J_num(3,2) = L3*(cosd(actual_angle(1,2))*cosd(actual_angle(1,3)) - sind(actual_angle(1,2))*sind(actual_angle(1,3))) + L2*cosd(actual_angle(1,2));
        J_num(1,3) = L3*(cosd(actual_angle(1,2))*sind(actual_angle(1,1))*sind(actual_angle(1,3)) + cosd(actual_angle(1,3))*sind(actual_angle(1,1))*sind(actual_angle(1,2)));
        J_num(2,3) = -L3*(cosd(actual_angle(1,1))*cosd(actual_angle(1,2))*sind(actual_angle(1,3)) + cosd(actual_angle(1,1))*cosd(actual_angle(1,3))*sind(actual_angle(1,2)));
        J_num(3,3) = L3*(cosd(actual_angle(1,2))*cosd(actual_angle(1,3)) - sind(actual_angle(1,2))*sind(actual_angle(1,3)));
                     
        % kolejne dane do nk
        A = J_num; 
        
     
       
        % najmniejsze kwadraty
        for i = 1:3
            for j = 1:3
                J_num(i,j) = 0;
                for z = 1:3
                    J_num(i,j) = J_num(i,j) + ( A(z,j) * weight(z,1) *  A(z,i));
                end
            end
        end
        
        
        
        %regularyzcja na diagonali

        for  i = 1:3 
            J_num(i,i) = J_num(i,i) + R(i,1);
        end
        
        
        
        %odwrocenie macierzy
        J_inv = inv(J_num);
        

        %ostatnie przemnozenie
        for i = 1:3 
            for j = 1:3
               J_num(i,j) = 0;
                for z = 1:3 
                    J_num(i,j) = J_num(i,j) + (J_inv(i,z)* A(j,z)  * weight(j,1));   
                end
            end
        end

        J_inv = J_num;


        %rozniczka polozen
        difer_x = poz_konc(1,1) - x_pos;
        difer_y = poz_konc(2,1) - y_pos;
        difer_z = poz_konc(3,1) - z_pos;


        %wektor rozniczki punktu
        diff_point = [difer_x;difer_y;difer_z];

        %wektor rozniczki kata
        diff_angle = [0;0;0];
        prev_diff_angle = [0;0;0];

        %obliczanie kroku przesuniecia kata
        for i = 1:3
           for j = 1:3  
               diff_angle(i,1) = diff_angle(i,1) + J_inv(i,j)*diff_point(j,1);
           end    
        end


        % dobranie aktualnego kata
        next_angle = [0;0;0];

        for i = 1:3
                next_angle(i,1) = actual_angle(1,i) + diff_angle(i,1);

                if (next_angle(i,1) > angleupperbound(i,1)) 
                    next_angle(i,1) = angleupperbound(i,1);
                end
                if (next_angle(i,1) < anglelowerbound(i,1)) 
                   next_angle(i,1) = anglelowerbound(i,1);
                end

                actual_angle(1,i) = next_angle(i,1);
                disp(actual_angle(1,i)); 

                serwo_theta_1 =  next_angle(1,1) + 90;
                serwo_theta_2 =  next_angle(2,1) + 90;
                serwo_theta_3 =  next_angle(3,1) + 180;
        end
        k = k + 1;
end





% petla ucinajaca niepotrzebna wielkosc z macierzy punktow
for i = k:10000
    punkty(k,:) = [];
    punkty_B(k,:) = [];
    punkty_C(k,:) = [];
end

figure(1);
hold on;
plot3(punkty(:,1),punkty(:,2),punkty(:,3),'-o');
plot3(punkty(1,1),punkty(1,2),punkty(1,3), 'mo');
plot3(0,0,0, 'mo');
plot3(poz_konc(1,1),poz_konc(2,1),poz_konc(3,1), 'mo');
for i = 1:k-1
    plot3(punkty_B(i,1),punkty_B(i,2),punkty_B(i,3), 'ro');
    plot3(punkty_C(i,1),punkty_C(i,2),punkty_C(i,3), 'go');
    line([0 punkty_C(i,1)],[0 punkty_C(i,2)],[0 punkty_C(i,3)],'Color','green','LineStyle','-');
    line([punkty_C(i,1) punkty_B(i,1)],[punkty_C(i,2) punkty_B(i,2)],[punkty_C(i,3) punkty_B(i,3)],'Color','red','LineStyle','-');
    line([punkty_B(i,1) punkty(i,1)],[punkty_B(i,2) punkty(i,2)],[punkty_B(i,3) punkty(i,3)],'Color','blue','LineStyle','-');
end    
xlabel('os X');
ylabel('os Y');
zlabel('os Z');
grid on;
hold off;


