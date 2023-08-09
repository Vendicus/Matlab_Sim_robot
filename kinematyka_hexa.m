%% program do kinematyki hexapoda na prace inzynierska.
% Aksamit Michał
%
clear;
clc;
close all;
%% wprowadzanie danych

theta_1 = 0; %48 ;
theta_2 = 0; %40 ;
theta_3 =  -90; %-123;

L1 = 51;
L2 = 70;
L3 = 207;
 
syms L_1
syms L_2
syms L_3
syms theta1
syms theta2
syms theta3

poz_konc = [-122;85.5;-207];
%% symboliczne

Rot_z0_sym = sym('O%d%d', [4 4]);

for i = 1:4
   for j = 1:4
       Rot_z0_sym(i,j) = 0;
   end    
end

Rot_z0_sym(1,1) = cos(theta1);
Rot_z0_sym(1,2) = -sin(theta1);
Rot_z0_sym(2,1) = sin(theta1);
Rot_z0_sym(2,2) = cos(theta1);
Rot_z0_sym(3,3) = 1;
Rot_z0_sym(4,4) = 1;


Rot_x0_sym = sym('O%d%d', [4 4]);

for i = 1:4
   for j = 1:4
       Rot_x0_sym(i,j) = 0;
   end    
end

Rot_x0_sym(1,1) = 1;
Rot_x0_sym(2,2) = cosd(-90);
Rot_x0_sym(3,2) = sind(-90);
Rot_x0_sym(2,3) = -sind(-90);
Rot_x0_sym(3,3) = cosd(-90);
Rot_x0_sym(4,4) = 1;


A0_sym = Rot_z0_sym*Rot_x0_sym;
A0_sym;


Trans_z1_sym =  sym('O%d%d', [4 4]);

for i = 1:4
   for j = 1:4
       Trans_z1_sym(i,j) = 0;
   end    
end

Trans_z1_sym(1,1) = 1;
Trans_z1_sym(2,2) = 1;
Trans_z1_sym(3,3) = 1;
Trans_z1_sym(4,4) = 1;

Trans_z1_sym(3,4) = L_1;


Rot_x1_sym = sym('O%d%d', [4 4]);

for i = 1:4
   for j = 1:4
       Rot_x1_sym(i,j) = 0;
   end    
end

Rot_x1_sym(1,1) = 1;
Rot_x1_sym(2,2) = cos(theta2);
Rot_x1_sym(3,2) = sin(theta2);
Rot_x1_sym(2,3) = -sin(theta2);
Rot_x1_sym(3,3) = cos(theta2);
Rot_x1_sym(4,4) = 1;

A1_sym = Trans_z1_sym * Rot_x1_sym;
A1_sym;


Trans_z2_sym =  sym('O%d%d', [4 4]);

for i = 1:4
   for j = 1:4
       Trans_z2_sym(i,j) = 0;
   end    
end

Trans_z2_sym(1,1) = 1;
Trans_z2_sym(2,2) = 1;
Trans_z2_sym(3,3) = 1;
Trans_z2_sym(4,4) = 1;

Trans_z2_sym(3,4) = L_2;

Rot_x2_sym = sym('O%d%d', [4 4]);

for i = 1:4
   for j = 1:4
       Rot_x2_sym(i,j) = 0;
   end    
end

Rot_x2_sym(1,1) = 1;
Rot_x2_sym(2,2) = cos(theta3);
Rot_x2_sym(3,2) = sin(theta3);
Rot_x2_sym(2,3) = -sin(theta3);
Rot_x2_sym(3,3) = cos(theta3);
Rot_x2_sym(4,4) = 1;

A2_sym = Trans_z2_sym * Rot_x2_sym;
A2_sym;

Trans_z3_sym =  sym('O%d%d', [4 4]);

for i = 1:4
   for j = 1:4
       Trans_z3_sym(i,j) = 0;
   end    
end

Trans_z3_sym(1,1) = 1;
Trans_z3_sym(2,2) = 1;
Trans_z3_sym(3,3) = 1;
Trans_z3_sym(4,4) = 1;

Trans_z3_sym(3,4) = L_3;

A3_sym = Trans_z3_sym;
A3_sym;

M_sym_B = A0_sym * A1_sym * A2_sym;
M_sym_C = A0_sym * A1_sym;
T_sym = A0_sym * A1_sym * A2_sym * A3_sym;

%% macierz A0

Rot_z0 = zeros(4);
Rot_x0 = zeros(4);

%przypisanie wartości do macierzy rotacji dla osi Z0
Rot_z0(1,1) = cosd(theta_1);
Rot_z0(1,2) = -sind(theta_1);
Rot_z0(2,1) = sind(theta_1);
Rot_z0(2,2) = cosd(theta_1);
Rot_z0(3,3) = 1;
Rot_z0(4,4) = 1;

%przypisanie wartości do macierzy rotacji dla osi X0
Rot_x0(1,1) = 1;
Rot_x0(2,2) = cosd(-90);
Rot_x0(3,2) = sind(-90);
Rot_x0(2,3) = -sind(-90);
Rot_x0(3,3) = cosd(-90);
Rot_x0(4,4) = 1;

%macierz przekształcenia A0 po przemnozeniu
A0 = Rot_z0 * Rot_x0;


%% macierz A1

Trans_z1 = zeros(4);

Trans_z1(1,1) = 1;
Trans_z1(2,2) = 1;
Trans_z1(3,3) = 1;
Trans_z1(4,4) = 1;

Trans_z1(3,4) = L1;

Rot_x1 = zeros(4);

Rot_x1(1,1) = 1;
Rot_x1(2,2) = cosd(theta_2);
Rot_x1(3,2) = sind(theta_2);
Rot_x1(2,3) = -sind(theta_2);
Rot_x1(3,3) = cosd(theta_2);
Rot_x1(4,4) = 1;


A1 = Trans_z1*Rot_x1;


%% macierz A2

Trans_z2 = zeros(4);

Trans_z2(1,1) = 1;
Trans_z2(2,2) = 1;
Trans_z2(3,3) = 1;
Trans_z2(4,4) = 1;

Trans_z2(3,4) = L2;

Rot_x2 = zeros(4);

Rot_x2(1,1) = 1;
Rot_x2(2,2) = cosd(theta_3);
Rot_x2(3,2) = sind(theta_3);
Rot_x2(2,3) = -sind(theta_3);
Rot_x2(3,3) = cosd(theta_3);
Rot_x2(4,4) = 1;

A2 = Trans_z2*Rot_x2;


%% macierz A3
Trans_z3 = zeros(4);

Trans_z3(1,1) = 1;
Trans_z3(2,2) = 1;
Trans_z3(3,3) = 1;
Trans_z3(4,4) = 1;

Trans_z3(3,4) = L3;

A3 = Trans_z3;


%% Macierz transformacji układu 

T = A0 * A1 * A2 * A3;
T


%% wykres 3D nogi 

%punkt C
M_C = A0*A1;
p_C = [M_C(1,4) M_C(2,4) M_C(3,4)];

%punkt B 
M_B = A0* A1* A2;
p_B = [M_B(1,4) M_B(2,4) M_B(3,4)];

%punkt A - efektor 
p_A = [T(1,4) T(2,4) T(3,4)];


figure(1);
plot3(M_C(1,4),M_C(2,4),M_C(3,4),'mo');
hold on;
plot3(M_B(1,4),M_B(2,4),M_B(3,4), 'mo');
plot3(T(1,4), T(2,4), T(3,4), 'mo');
plot3(0,0,0,'mo');
line([0 p_C(1,1)],[0 p_C(1,2)],[0 p_C(1,3)],'Color','blue','LineStyle','-');
line([p_C(1,1) p_B(1,1)],[p_C(1,2) p_B(1,2)],[p_C(1,3) p_B(1,3)],'Color','green','LineStyle','-');
line([p_B(1,1) p_A(1,1)],[p_B(1,2) p_A(1,2)],[p_B(1,3) p_A(1,3)],'Color','red','LineStyle','-');
grid on;
xlabel('x');
ylabel('y');
zlabel('z');

text(M_C(1,4),M_C(2,4),M_C(3,4)+20,'C', 'color', 'blue');
text(M_B(1,4),M_B(2,4),M_B(3,4)+20,'B', 'color', 'green');
text(p_A(1,1),p_A(1,2),p_A(1,3)-20,'A', 'color', 'red');
text(0,0,20,'D', 'color', 'black');
hold off;


%% jakobiany 

J = jacobian([T_sym(1,4),T_sym(2,4),T_sym(3,4)],[theta1, theta2, theta3]);

% odwrotność jakobiana
J_inv = inv(J);
%transpozycja jakobiana
J_trans = transpose(J);
%tworzenie jakobiana do obliczeń numerycznych
J_num = zeros(3);

J_num(1,1) = L3*(cosd(theta_1)*sind(theta_2)*sind(theta_3) - cosd(theta_1)*cosd(theta_2)*cosd(theta_3)) - L1*cosd(theta_1) - L2*cosd(theta_1)*cosd(theta_2);
J_num(2,1) = L3*(sind(theta_1)*sind(theta_2)*sind(theta_3) - cosd(theta_2)*cosd(theta_3)*sind(theta_1)) - L1*sind(theta_1) - L2*cosd(theta_2)*sind(theta_1);
J_num(3,1) = 0;
J_num(1,2) = L3*(cosd(theta_2)*sind(theta_1)*sind(theta_3) + cosd(theta_3)*sind(theta_1)*sind(theta_2)) + L2*sind(theta_1)*sind(theta_2);
J_num(2,2) = -L3*(cosd(theta_1)*cosd(theta_2)*sind(theta_3) + cosd(theta_1)*cosd(theta_3)*sind(theta_2)) - L2*cosd(theta_1)*sind(theta_2);
J_num(3,2) = L3*(cosd(theta_2)*cosd(theta_3) - sind(theta_2)*sind(theta_3)) + L2*cosd(theta_2);
J_num(1,3) = L3*(cosd(theta_2)*sind(theta_1)*sind(theta_3) + cosd(theta_3)*sind(theta_1)*sind(theta_2));
J_num(2,3) = -L3*(cosd(theta_1)*cosd(theta_2)*sind(theta_3) + cosd(theta_1)*cosd(theta_3)*sind(theta_2));
J_num(3,3) = L3*(cosd(theta_2)*cosd(theta_3) - sind(theta_2)*sind(theta_3));

% odwrocenie macierzy Jakobianu num
J_num_inv = inv(J_num)

%transpozycja macierzy Jakobianu num
J_num_trans = transpose(J_num)

%% kroki do kątów

%roznica miedzy pozycjami startowymi, koncowymi
difer_x = sqrt((T(1,4) - poz_konc(1,1))^2);
difer_y = sqrt((T(2,4) - poz_konc(2,1))^2);
difer_z = sqrt((T(3,4) - poz_konc(3,1))^2);

diff_point = [difer_x;difer_y;difer_z];

diff_angle = [0;0;0];
prev_diff_angle = [0;0;0];

%obliczanie kroku przesuniecia kata
for i = 1:3
   for j = 1:3
       diff_angle(i,1) = prev_diff_angle(i,1) + J_num_inv(i,j)*diff_point(j,1);
       prev_diff_angle(i,1) = diff_angle(i,1);
   end    
end


% dobranie aktualnego kata
actual_angle = [theta_1;theta_2;theta_3];
next_angle = [0;0;0];

for i = 1:3
		next_angle(i,1) = actual_angle(i,1) + diff_angle(i,1);
%{
        if (next_angle(i,1) > angleupperbound[ii+1]) 
            prevangle[ii+1] = angleupperbound[ii+1];
        end
		if (next_angle(i,1) < anglelowerbound[ii+1]) 
            prevangle[ii+1] = anglelowerbound[ii+1];
        end
        
%}    
        
        actual_angle(i,1) = next_angle(i,1);
end

