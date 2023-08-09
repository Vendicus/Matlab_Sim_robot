%% program do fizyki hexapoda na prace inzynierska.
% Aksamit Michał
%

clc; 
close all; 
clear variables;
%% dane:

% podaj dane długości członow hexapoda ( mm )

L1 = 207;   % golen
L2 = 70;    % udo
L3 = 51;    % biodro
LK = 80;    % odleglosc od srodka ciezkosci do biodra

% Podaj masy danych modulow (g) :

M_bateria = 451;
M_elektronika = 87;
M_szkielet = 850;

M_serwo = 43;

M_czlon_1 = 49;
M_czlon_2 = 25;
M_czlon_3 = 13;

%podaj kat zgiecia do osi x dla uda ( stopnie )
theta_2 = 10;

%podaj kat zgiecia do osi x dla goleni ( stopnie )
theta_1 = 60;

%podaj kat rozwarcia nog od osi x dla lewej strony robota ( stopnie )
kat_rozw_pocz = -20;
kat_rozw_konc = 20;
krok = 1;

kat_rozw = kat_rozw_pocz:krok:kat_rozw_konc;

%kat rozwarcia nog srodkowych od osi x ( stopnie )
kat_prawa = 0 ;%kat_rozw_pocz:krok:kat_rozw_konc;

%w celu symulacji geometrii należy przyjac jeden kat rozw bez przedzialu
%UWAGA: GDY SYMULACJA AKTYWNA, nie mozna skorzystac z wykresow momentow i
%predkosci katowych, natomiast gdy symulacja nie aktywna, nie mozna zobaczyc geometrii.

kat_rozw_sym = 30;
kat_prawa_sym = 0;

symulacja = 0; % czy symulowac? (1 - tak, 0 - nie)

%czy robot w fazie stania czy kroczenia? (1 - kroczenie, 0 - stanie)
kroczenie = 1;

%predkosc serwo (deg/s) bez obciazenia dla 4,8 V
w = 375;
%moment serwa [kgcm]:
Mom_serwo = 5.8;

%dane dla zrodlo zasilania:
pojemnosc_ba = 8;   %AH
na_ba = 11.1;  %V
C_ba = 30; %rozladowanie

spr_rp = 80/100 ; %sprawnosc przetwornicy dla danej wartosci nspiecia Rpi.
spr_sters = 95/100 ; %sprawnosc przetwornicy dla danej wartosci napiec na sterownik serw
spr_bateria = 90/100; %sprawnosc baterii LI-PO

%-- dla serwa powerhd 6001 HB 
%Pobór prądu w stanie spoczynku: 4 mA
%Pobór prądu bez obciążenia: 250 mA
%Pobór prądu maksymalny: 1300 mA

%podaj przyspieszenie ziemskie ( m/s^2 )
g = 9.8105;


%% konwersja jednostek do SI, liczenie sil ciezkosci czlonow

% konwersja do metr
LK = LK * 0.001;
L1 = L1 * 0.001;
L2 = L2 * 0.001;
L3 = L3 * 0.001;

% konwersja do kg
M_bateria = M_bateria * 0.001;
M_elektronika = M_elektronika * 0.001;
M_szkielet = M_szkielet * 0.001;
M_serwo = M_serwo * 0.001;
M_czlon_1 = M_czlon_1 * 0.001;
M_czlon_2 = M_czlon_2 * 0.001;
M_czlon_3 = M_czlon_3 * 0.001;

%obliczanie masy nogi
M_nogi = (M_czlon_1 + M_czlon_2 + M_czlon_3) + (3 * M_serwo);
fprintf(" Masa nogi robota wynosi : %.4f  kg \n",M_nogi);

%obliczanie masy korpusu
M_korpus = (M_bateria + M_elektronika + M_szkielet) + (3 * M_nogi);
fprintf(" Masa korpusu robota (z 3 podniesionymi nogami) wynosi : %.4f  kg \n",M_korpus);

%obliczanie masy calkowitej
M_calkowita = M_korpus + 3 * M_nogi;
fprintf(" Masa calkowita robota wynosi : %.4f  kg \n\n", M_calkowita);



%------------ liczenie sil ciezkosci --------------

Q_korpus = M_korpus * g;
Q_noga = M_nogi * g;
Q_1 = (M_czlon_1 + M_serwo) * g;
Q_2 = (M_czlon_2 + M_serwo) * g;
Q_3 = (M_czlon_3 + M_serwo) * g;


%% obliczanie czasu pracy na baterii;

moc_baterii = pojemnosc_ba * na_ba;

%pobor mocy Rpi

moc_Rpi = 5 * 3 * (1/spr_rp);

%zapotrzebowanie na moc na jeden servobonnet
%zakladam 3 serwa - duże (30% kgcm dostepnych na serwo) obciazenia, 6 serw - obciazenia mniejsze ( 25% ) 
ampery_serwo = 3*0.5 + 6*0.4;
moc_serwo_bonnet = ampery_serwo *5*(1/spr_sters);

%zapotrzebowanie na moc pelnego ukladu

moc_ukl = 2*moc_serwo_bonnet + moc_Rpi;
ampery_ukl = 3 * (1/spr_rp) + ampery_serwo*(1/spr_sters);

%czas pracy
czas_pracy_teor = moc_baterii/moc_ukl;

% zalozenie sprawnosci baterii i odlaczenia przy naładowaniu wynoszącym
% tylko 20%.
czas_pracy_rzecz = czas_pracy_teor * spr_bateria *(1-(20/100));

fprintf(" Rzeczywisty czas pracy bateri wynosi : %.4f h, a teoretyczny wynosi : %.4f h \n\n",czas_pracy_rzecz, czas_pracy_teor);


%% obliczanie sil reakcji z twierdzen rownowagi momentow i rownowagi sil - STATYKA

% reakcja R2 z  rownania momentu statycznego wzg. punktu A. Reakcja ta opisuje sile reakcji R2 na osi y na stopie
% lewej.
R2_y = ( Q_1*L1*cosd(theta_1)*cosd(kat_prawa) + Q_2*( L1*cosd(theta_1) + L2*cosd(theta_2) )*cosd(kat_prawa) + Q_3 * ( L3 + L2*cosd(theta_2) + L1*cosd(theta_1) )*cosd(kat_prawa) + Q_korpus * ( LK + ( L3 + L2*cosd(theta_2) + L1*cosd(theta_1) )*cosd(kat_prawa) ) + 2*Q_3*( 1.5*LK + ( L3 + L2*cosd(theta_2) + L1*cosd(theta_1) )*cosd(kat_prawa) ) + 2*Q_2*( 1.5*LK + L3*cosd(kat_rozw) + ( L3 + L2*cosd(theta_2) + L1*cosd(theta_1) )*cosd(kat_prawa) ) + 2*Q_1*( 1.5*LK + L3*cosd(kat_rozw) + L2*cosd(theta_2)*cosd(kat_rozw) +  ( L3 + L2*cosd(theta_2) + L1*cosd(theta_1) )*cosd(kat_prawa) ) ) / ( 2*( 1.5*LK + L3*cosd(kat_rozw) + L3*cosd(kat_prawa) + L2*cosd(theta_2)*cosd(kat_rozw) + L2*cosd(theta_2)*cosd(kat_prawa) + L1*cosd(theta_1)*cosd(kat_rozw) + L1*cosd(theta_1)*cosd(kat_prawa) ) );
R2 = R2_y / sind(theta_1);
R2_x = R2 * cosd(theta_1);
fprintf(" Reakcja R2 robota wynosi : %.4f  N \n",R2);

% obliczanie reakcji R2 z rownania rownowagi sil dla robota na osi y
R1_y = Q_korpus + (3 * Q_noga) - 2*R2_y;
R1 = R1_y / sind(theta_1);
R1_x = R1*cosd(theta_1);
fprintf(" Reakcja R1 robota wynosi : %.4f  N \n\n\n",R1);



% ------------------ MOMENTY na serwo w okreslonych punktach przegubow członow --------------------------

% obliczanie sumy momentow wzg punktu B
M_B = -( R1_y * L1 * cosd(theta_1) * cosd(kat_prawa) - (R1_x * L1*sind(theta_1)*cosd(kat_prawa)) + Q_2*L2*cosd(theta_2)*cosd(kat_prawa) + Q_3*(L3+L2*cosd(theta_2))*cosd(kat_prawa) + Q_korpus*(LK + ((L3 + L2*cosd(theta_2))*cosd(kat_prawa))) + 2*Q_3*( (1.5 * LK) + L3*cosd(kat_prawa) + L2*cosd(theta_2)*cosd(kat_prawa)) + 2*Q_2*(1.5*LK + (cosd(kat_rozw)*L3 + L3*cosd(kat_prawa)) + L2*cosd(theta_2)*cosd(kat_prawa)) + 2*Q_1*( 1.5*LK + L3*cosd(kat_prawa) + L3*cosd(kat_rozw) + L2*cosd(theta_2)*cosd(kat_prawa) + L2*cosd(theta_2)*cosd(kat_rozw) )  - 2*R2_y*(1.5*LK + cosd(kat_rozw)*L3 + L3*cosd(kat_prawa) + L2*cosd(theta_2)*cosd(kat_prawa) + L2*cosd(theta_2)*cosd(kat_rozw) + L1*cosd(theta_1)*cosd(kat_rozw)) + 2*R2_x*L1*sind(theta_1)*cosd(kat_rozw));
M_B = M_B * 10.1972; % konwersja na Kgcm z Nm
fprintf(" Moment wzg punktu B laczenia robota wynosi : %.4f  kgcm \n",M_B);
fprintf("\n");

% obliczanie sumy momentow wzg punktu C
M_C = (Q_1 * L2*cosd(theta_2)*cosd(kat_prawa)) - (R1_y * cosd(kat_prawa) * (L2*cosd(theta_2) + L1*cosd(theta_1))) + R1_x * cosd(kat_prawa)*( L1*sind(theta_1) - L2*sind(theta_2)) - Q_3*L3*cosd(kat_prawa) - Q_korpus*( LK + L3*cosd(kat_prawa)) - 2*Q_3*( 1.5*LK + L3*cosd(kat_prawa) ) -  2*Q_2*(1.5*LK + (cosd( kat_rozw )*L3 + L3*cosd(kat_prawa))) -  2*Q_1*(1.5*LK + (L3*cosd(kat_prawa) + cosd(kat_rozw)*L3) + cosd(kat_rozw)*L2*cosd(theta_2)) + 2*R2_y*(cosd(kat_rozw)*L1*cosd(theta_1) + cosd(kat_rozw)*L2*cosd(theta_2) + 1.5*LK + (cosd(kat_rozw)*L3 + L3*cosd(kat_prawa))) - 2*R2_x*cosd(kat_rozw)*( L1*sind(theta_1) - L2*sind(theta_2));
M_C = M_C * 10.1972;
fprintf(" Moment wzg punktu C laczenia robota wynosi : %.4f  kgcm \n",M_C);
fprintf("\n");

% obliczanie sumy momentow wzg punktu D
M_D = -(R1_y * (L3 + L2*cosd(theta_2) + L1*cosd(theta_1))*cosd(kat_prawa)) + R1_x * cosd(kat_prawa) * ( L1*sind(theta_1) - L2*sind(theta_2)) + Q_1*cosd(kat_prawa)*(L3 + L2*cosd(theta_2)) + Q_2*L3*cosd(kat_prawa) -(Q_korpus*LK) -2*Q_3*(1.5*LK) - 2*Q_2*(1.5*LK + (cosd(kat_rozw)*L3)) - 2*Q_1*((cosd(kat_rozw)*(L2*cosd(theta_2)+L3))+1.5*LK) + 2*R2_y*( cosd(kat_rozw)*L1*cosd(theta_1) + cosd(kat_rozw)*L2*cosd(theta_2) + cosd(kat_rozw)*L3 + 1.5*LK) - 2*R2_x*cosd(kat_rozw)*( L1*sind(theta_1) - L2*sind(theta_2));  
M_D = M_D * 10.1972;
fprintf(" Moment wzg punktu D laczenia robota wynosi : %.4f  kgcm \n",M_D);
fprintf("\n\n");

%% wykresy momentow i maksymalne/minimalne momenty

% nalezy podac zakres dla katow rozwarcia w danych np. -20:1:20 (od -20 stopni do 20 stopni z krokiem 1).


if( symulacja == 0)

    
    
if( sqrt(min(M_B)^2) >= sqrt(max(M_B)^2) )
    [maks_B, maks_XB] = min(M_B);
else
    [maks_B, maks_XB] = max(M_B);
end    

if( sqrt(min(M_C)^2) >= sqrt(max(M_C)^2) )
    [maks_C, maks_XC] = min(M_C);
else
    [maks_C, maks_XC] = max(M_C);
end    

if( sqrt(min(M_D)^2) >= sqrt(max(M_D)^2) )
    [maks_D, maks_XD] = min(M_D);
else
    [maks_D, maks_XD] = max(M_D);
end    

figure();
subplot(2,2,[1,2]);
plot( kat_rozw, M_B,'LineWidth',1.5, 'Color', 'r');
title(' Wykres momentów w punkcie B. ', 'FontSize', 15);
xlabel('kąt rozwarcia nóg ( kąt_{rozw}  i  kąt_{prawa} ) [\circ]','FontSize', 12, 'Color', 'black');
ylabel(' Moment [Kgcm] ', 'FontSize', 12, 'Color', 'black');
hold on;
grid on;
plot((maks_XB + (kat_rozw_pocz - 1)*krok), maks_B, 'o', 'MarkerSize',10, 'MarkerEdgeColor','y', 'MarkerFaceColor','white');
text((maks_XB + (kat_rozw_pocz - 1)*krok), maks_B+0.01,[' Największy moment dla punktu B wynosi : ' num2str(maks_B) ' Kgcm']);
hold off;

subplot(2,2,3);
plot( kat_rozw, M_C,'LineWidth',1.5, 'Color', 'g');
title(' Wykres momentów w punkcie C. ', 'FontSize', 15);
xlabel('kąt rozwarcia nóg ( kąt_{rozw}  i  kąt_{prawa} ) [\circ]','FontSize', 12, 'Color', 'black');
ylabel(' Moment [Kgcm] ', 'FontSize', 12, 'Color', 'black');
hold on;
grid on;
plot((maks_XC + (kat_rozw_pocz - 1)*krok), maks_C, 'o', 'MarkerSize',10, 'MarkerEdgeColor','b', 'MarkerFaceColor','white');
text((maks_XC + (kat_rozw_pocz - 1)*krok), maks_C+0.01,[' Największy moment dla punktu C wynosi : ' num2str(maks_C) ' Kgcm']);
hold off;

subplot(2,2,4);
plot(kat_rozw, M_D,'LineWidth',1.5, 'Color', 'b');
title(' Wykres momentów w punkcie D. ', 'FontSize', 15);
xlabel('kąt rozwarcia nóg ( kąt_{rozw}  i  kąt_{prawa} ) [\circ]','FontSize', 12, 'Color', 'black');
ylabel(' Moment [Kgcm] ', 'FontSize', 12, 'Color', 'black');
hold on;
grid on;
plot((maks_XD + (kat_rozw_pocz - 1)*krok), maks_D, 'o', 'MarkerSize',10, 'MarkerEdgeColor','cyan', 'MarkerFaceColor','white');
text((maks_XD + (kat_rozw_pocz - 1)*krok), maks_D+0.01,[' Największy moment dla punktu D wynosi : ' num2str(maks_D) ' Kgcm']);
hold off;


%% kinematyka przykładowej nogi prawej.

% ---------------------------- obliczanie predkosci katowej i przyspieszenia wzgledem punktu D.

w_rzecz_D = ((Mom_serwo - abs(M_D))/Mom_serwo) * w;
t_w_D = zeros(1,length(w_rzecz_D)); % tworzenie tablicy czasu do obliczen zaleznosci predkosci katowej od czasu

for i = 1:1:(length(w_rzecz_D)-1)     %nadawanie wartosci do tablicy czasu
    if(i-1 == 0)
        t_w_D(i+1) = ((2 * krok) / ( w_rzecz_D(i) + w_rzecz_D(i+1) ));  %obliczanie czasu uzyskania danego kata uwzgledniajac predkosc katowa
    else
        t_w_D(i+1) = ((2 * krok) / ( w_rzecz_D(i) + w_rzecz_D(i+1) )) + t_w_D(i);
    end
end

e_D = zeros(1,length(w_rzecz_D));

for i = 1:1:(length(w_rzecz_D)-1)     %nadawanie wartosci do tablicy czasu    
            e_D(i+1) = ( w_rzecz_D(i+1) - w_rzecz_D(i) )/( t_w_D(i+1) - t_w_D(i) ); %obliczanie przyspieszenia katowego
end

figure();
subplot(2,2,1);
plot( t_w_D, w_rzecz_D, 'LineWidth', 1.5, 'Color', 'r');
title(' Prędkość kątowa \omega  wzgl. punktu D w zależnośći od czasu', 'FontSize', 15);
xlabel(' czas [s] ','FontSize', 12, 'Color', 'black');
ylabel(' prędkość kątowa \omega [\circ/s] ', 'FontSize', 12, 'Color', 'black');
grid on;

subplot(2,2,2);
plot( t_w_D, kat_prawa, 'LineWidth', 1.5, 'Color', 'g');
title(' kąt \theta_{2} w zależnośći od czasu wzgl. punktu D ', 'FontSize', 15);
xlabel(' czas [s] ','FontSize', 12, 'Color', 'black');
ylabel(' kąt [\circ] ', 'FontSize', 12, 'Color', 'black');
grid on;


subplot(2,2,3);
plot( kat_rozw, w_rzecz_D, 'LineWidth', 1.5, 'Color', 'r');
title(' Prędkość \omega wzg. punktu D w zależnośći od kąta_{ rozw} ', 'FontSize', 15);
xlabel(' kąt rozwarcia [\circ] ','FontSize', 12, 'Color', 'black');
ylabel(' prędkość kątowa \omega [\circ/s] ', 'FontSize', 12, 'Color', 'black');
grid on;


subplot(2,2,4);
plot( t_w_D, e_D, 'LineWidth', 1.5, 'Color', 'b');
title(' przyśpieszenie kątowe \epsilon wzg. punktu D w zależnośći od czasu ', 'FontSize', 15);
xlabel(' czas [s] ','FontSize', 12, 'Color', 'black');
ylabel(' przyśpieszenie kątowe \epsilon [\circ/s^{2}] ', 'FontSize', 12, 'Color', 'black');
grid on;




% ----------------------------------- predkosc katowa i przyspieszenie katowe wzg. punktu C

w_rzecz_C = ((Mom_serwo - abs(M_C))/Mom_serwo) * w;
t_w_C = zeros(1,length(w_rzecz_C)); % tworzenie tablicy czasu do obliczen zaleznosci predkosci katowej od czasu

for i = 1:1:(length(w_rzecz_C)-1)     %nadawanie wartosci do tablicy czasu
    if(i-1 == 0)
        t_w_C(i+1) = ((2 * krok) / ( w_rzecz_C(i) + w_rzecz_C(i+1) ));  %obliczanie czasu uzyskania danego kata uwzgledniajac predkosc katowa
    else
        t_w_C(i+1) = ((2 * krok) / ( w_rzecz_C(i) + w_rzecz_C(i+1) )) + t_w_C(i);
    end
end

e_C = zeros(1,length(w_rzecz_C));

for i = 1:1:(length(w_rzecz_C)-1)     %nadawanie wartosci do tablicy czasu    
            e_C(i+1) = ( w_rzecz_C(i+1) - w_rzecz_C(i) )/( t_w_C(i+1) - t_w_C(i) ); %obliczanie przyspieszenia katowego
end

figure();
subplot(2,2,1);
plot( t_w_C, w_rzecz_C, 'LineWidth', 1.5, 'Color', 'r');
title(' Prędkość kątowa \omega  wzgl. punktu C w zależnośći od czasu', 'FontSize', 15);
xlabel(' czas [s] ','FontSize', 12, 'Color', 'black');
ylabel(' prędkość kątowa \omega [\circ/s] ', 'FontSize', 12, 'Color', 'black');
grid on;

subplot(2,2,2);
plot( t_w_C, kat_prawa, 'LineWidth', 1.5, 'Color', 'g');
title(' kąt \theta_{2} w zależnośći od czasu wzgl. punktu C ', 'FontSize', 15);
xlabel(' czas [s] ','FontSize', 12, 'Color', 'black');
ylabel(' kąt [\circ] ', 'FontSize', 12, 'Color', 'black');
grid on;


subplot(2,2,3);
plot( kat_rozw, w_rzecz_C, 'LineWidth', 1.5, 'Color', 'r');
title(' Prędkość \omega wzg. punktu C w zależnośći od kąta_{ rozw} ', 'FontSize', 15);
xlabel(' kąt rozwarcia [\circ] ','FontSize', 12, 'Color', 'black');
ylabel(' prędkość kątowa \omega [\circ/s] ', 'FontSize', 12, 'Color', 'black');
grid on;


subplot(2,2,4);
plot( t_w_C, e_C, 'LineWidth', 1.5, 'Color', 'b');
title(' przyśpieszenie kątowe \epsilon wzg. punktu C w zależnośći od czasu ', 'FontSize', 15);
xlabel(' czas [s] ','FontSize', 12, 'Color', 'black');
ylabel(' przyśpieszenie kątowe \epsilon [\circ/s^{2}] ', 'FontSize', 12, 'Color', 'black');
grid on;



%  ------------------------------------------ predkosc katowa i przyspieszenie katowe wzg. punktu B

w_rzecz_B = ((Mom_serwo - abs(M_B))/Mom_serwo) * w;
t_w_B = zeros(1,length(w_rzecz_B)); % tworzenie tablicy czasu do obliczen zaleznosci predkosci katowej od czasu

for i = 1:1:(length(w_rzecz_B)-1)     %nadawanie wartosci do tablicy czasu
    if(i-1 == 0)
        t_w_B(i+1) = ((2 * krok) / ( w_rzecz_B(i) + w_rzecz_B(i+1) ));  %obliczanie czasu uzyskania danego kata uwzgledniajac predkosc katowa
    else
        t_w_B(i+1) = ((2 * krok) / ( w_rzecz_B(i) + w_rzecz_B(i+1) )) + t_w_B(i);
    end
end

e_B = zeros(1,length(w_rzecz_B));

for i = 1:1:(length(w_rzecz_B)-1)     %nadawanie wartosci do tablicy czasu    
            e_B(i+1) = ( w_rzecz_B(i+1) - w_rzecz_B(i) )/( t_w_B(i+1) - t_w_B(i) ); %obliczanie przyspieszenia katowego
end

figure();
subplot(2,2,1);
plot( t_w_B, w_rzecz_B, 'LineWidth', 1.5, 'Color', 'r');
title(' Prędkość kątowa \omega  wzgl. punktu B w zależnośći od czasu', 'FontSize', 15);
xlabel(' czas [s] ','FontSize', 12, 'Color', 'black');
ylabel(' prędkość kątowa \omega [\circ/s] ', 'FontSize', 12, 'Color', 'black');
grid on;

subplot(2,2,2);
plot( t_w_B, kat_prawa, 'LineWidth', 1.5, 'Color', 'g');
title(' kąt \theta_{2} w zależnośći od czasu wzgl. punktu B ', 'FontSize', 15);
xlabel(' czas [s] ','FontSize', 12, 'Color', 'black');
ylabel(' kąt [\circ] ', 'FontSize', 12, 'Color', 'black');
grid on;


subplot(2,2,3);
plot( kat_rozw, w_rzecz_B, 'LineWidth', 1.5, 'Color', 'r');
title(' Prędkość \omega wzgl. pukntu B w zależnośći od kąta \theta_{2} ', 'FontSize', 15);
xlabel(' kąt rozwarcia goleni ( \theta_{2} ) [\circ] ','FontSize', 12, 'Color', 'black');
ylabel(' prędkość kątowa \omega [\circ/s] ', 'FontSize', 12, 'Color', 'black');
grid on;


subplot(2,2,4);
plot( t_w_B, e_B, 'LineWidth', 1.5, 'Color', 'b');
title(' przyśpieszenie kątowe \epsilon wzgl. punktu B w zależnośći od czasu ', 'FontSize', 15);
xlabel(' czas [s] ','FontSize', 12, 'Color', 'black');
ylabel(' przyśpieszenie kątowe \epsilon [\circ/s^{2}] ', 'FontSize', 12, 'Color', 'black');
grid on;

end



%% geometria
skala = 100;


if ( symulacja == 1)
    

kat_rozw = kat_rozw_sym;
kat_prawa = kat_prawa_sym;


%wspolrzedne punktu D
XD = LK *skala;
YD = 0;
D = [XD YD];

%wspolrzedne punktu C
XC = ( LK + L3*cosd(kat_prawa) ) *skala;
YC = 0;
C = [XC YC];

%wspolrzedne punktu B
XB = ( LK + (L3 + L2*cosd(theta_2))*cosd(kat_prawa)) *skala;
YB = (L2*sind(theta_2)) *skala;
B = [XB YB];

%wspolrzedne punktu A
XA = ( LK + (L3 + L2*cosd(theta_2) + L1*cosd(theta_1))*cosd(kat_prawa)) *skala;
YA = (L2*sind(theta_2) - L1*sind(theta_1)) *skala;
A = [XA YA];

%wspolrzedne punktu E
XE = -(LK * 0.5)*skala;
YE = 0;

%wspolrzedne punktu F
XF = -(((LK * 0.5) + L3*cosd(kat_rozw) )*skala);
YF = 0;

%wspolrzedne punktu G
XG = -(((LK * 0.5) + L3*cosd(kat_rozw) + L2*cosd(theta_2)*cosd(kat_rozw) )*skala);
YG = YB;

%wspolrzedne punktu H
XH = -(((LK * 0.5) + L3*cosd(kat_rozw) + L2*cosd(theta_2)*cosd(kat_rozw) +L1*cosd(theta_1)*cosd(kat_rozw) )*skala);
YH = YA;


limit_x = sqrt((XA * 1.35)^2);
limit_y = sqrt((YA * 2.25)^2);

%zdefiniowanie punktów przegubów robota i jego rzutu z przodu
line([XH XG XF XE XD XC XB XA],[YH YG YF YE YD YC YB YA], 'Color', 'b', 'LineWidth', 1.5);
title(' Rzut z przodu na robota z czlonami podpierajacymi ', 'FontSize', 15);
xlim([-limit_x limit_x]);
ylim([-limit_y limit_y]);
hold on;
plot(XH,YH,'o', 'MarkerSize',10, 'MarkerEdgeColor','r', 'MarkerFaceColor','white');
text(XH -(0.001*skala),YH - (0.02*skala),'H','Color','red','FontSize',12,  'HorizontalAlignment','right');
plot(XG,YG,'o', 'MarkerSize',10, 'MarkerEdgeColor','c', 'MarkerFaceColor','white');
text(XG -(0.001*skala),YG + (0.02*skala),'G','Color','cyan','FontSize',12,  'HorizontalAlignment','right');
plot(XF,YF,'o', 'MarkerSize',10, 'MarkerEdgeColor','g', 'MarkerFaceColor','white');
text(XF -(0.001*skala),YF + (0.02*skala),'F','Color','green','FontSize',12,  'HorizontalAlignment','right');
plot(XE,YE,'o', 'MarkerSize',10, 'MarkerEdgeColor','b', 'MarkerFaceColor','white');
text(XE -(0.001*skala),YE + (0.02*skala),'E','Color','blue','FontSize',12,  'HorizontalAlignment','right');
plot(XD,YD,'o', 'MarkerSize',10, 'MarkerEdgeColor','b', 'MarkerFaceColor','white');
text(XD +(0.001*skala),YD + (0.02*skala),'D','Color','blue','FontSize',12,  'HorizontalAlignment','left');
plot(XC,YC,'o', 'MarkerSize',10, 'MarkerEdgeColor','g', 'MarkerFaceColor','white');
text(XC +(0.001*skala),YC + (0.02*skala),'C','Color','green','FontSize',12,  'HorizontalAlignment','left');
plot(XB,YB,'o', 'MarkerSize',10, 'MarkerEdgeColor','c', 'MarkerFaceColor','white');
text(XB +(0.001*skala),YB + (0.02*skala),'B','Color','cyan','FontSize',12,  'HorizontalAlignment','left');
plot(XA,YA,'o', 'MarkerSize',10, 'MarkerEdgeColor','r', 'MarkerFaceColor','white');
text(XA +(0.001*skala),YA - (0.02*skala),'A','Color','red','FontSize',12, 'HorizontalAlignment','left');
plot(0,0,'square', 'MarkerSize',15, 'MarkerEdgeColor','black', 'MarkerFaceColor','white');
text(0,0 + (0.02*skala),'śr. ciężkości korpusu','Color','black','FontSize',13,  'HorizontalAlignment','center');
grid on;

%ziemia
yline(YA,':','LineWidth', 1.5);


% ---- katy ----

%definiuj kat theta-1
rozdz_th1 = ( 180 - theta_1):0.1:180;

%centrum kata
r = skala*0.05;

% def rownania
X_theta_1 = XA + r*cosd(rozdz_th1);
Y_theta_1 = YA + r*sind(rozdz_th1);

plot(X_theta_1,Y_theta_1, 'Color', 'black','LineStyle','-', 'LineWidth', 0.5);
text((XA/1.2),(YA+(skala*0.015)),['\it \theta_{1} = ' num2str(theta_1) '\circ'], 'Color', 'black', 'FontSize',13, 'HorizontalAlignment','left');


%def kat theta-2

rozdz_th2 = 0:0.1:theta_2;
r_2 = skala*0.05;
X_theta_2 = XC + r_2*cosd(rozdz_th2);
Y_theta_2 = YC + r_2*sind(rozdz_th2);
plot(X_theta_2,Y_theta_2, 'Color', 'black','LineStyle','-', 'LineWidth', 0.5);

if(theta_2 >= 30)
    text((XC*1.2),(YC+(skala*0.0035)),['\it \theta_{2} = ' num2str(theta_2) '\circ'], 'Color', 'black', 'FontSize',13, 'HorizontalAlignment','left');
elseif(theta_2 < 30 && theta_2 > 0)
    text((XC*1.2),(YC-(skala*0.002)),['\it \theta_{2} = ' num2str(theta_2) '\circ'], 'Color', 'black', 'FontSize',13, 'HorizontalAlignment','left');      
else
    text((XC*1.2),(YC-(skala*0.005)),['\it \theta_{2} = ' num2str(theta_2) '\circ'], 'Color', 'black', 'FontSize',13, 'HorizontalAlignment','left');   
end

line([ XC (X_theta_2(end) + skala*0.01)],[YC YC], 'Color', 'Black', 'LineStyle', '--', 'LineWidth', 0.5);


% ---------- dlugosci czlonow ----------


%L1
plot([XA XB],[(YB+(skala*0.1)) (YB+(skala*0.1))],'-','LineWidth',0.5,'Color','black');
plot([XA XB],[(YB+(skala*0.1)) (YB+(skala*0.1))],'|','LineWidth',0.5,'Color','black');
text(((XA+XB)/2),(YB+(skala*0.12)),['\it L_{1}cos(' num2str(theta_1) '\circ) = ' num2str(L1*cosd(theta_1)*cosd(kat_prawa)*skala) ' cm'], 'Color', 'black', 'FontSize',11, 'HorizontalAlignment','left');
plot([XA XA], [YA (YB+(skala*0.1))], '-','LineWidth',0.1,'Color','black');
plot([XB XB], [YB (YB+(skala*0.1))], '-','LineWidth',0.1,'Color','black');
%L2
plot([XB XC],[(YB+(skala*0.1)) (YB+(skala*0.1))],'-','LineWidth',0.5,'Color','black');
plot([XB XC],[(YB+(skala*0.1)) (YB+(skala*0.1))],'|','LineWidth',0.5,'Color','black');
text(((XB+XC)/2),(YB+(skala*0.12)),['\it L_{2}cos(' num2str(theta_2) '\circ) = ' num2str(L2*cosd(theta_2)*cosd(kat_prawa)*skala) ' cm'], 'Color', 'black', 'FontSize',11, 'HorizontalAlignment','center');
plot([XC XC], [YC (YB+(skala*0.1))], '-','LineWidth',0.1,'Color','black');
%L3
plot([XC XD],[(YB+(skala*0.1)) (YB+(skala*0.1))],'-','LineWidth',0.5,'Color','black');
plot([XC XD],[(YB+(skala*0.1)) (YB+(skala*0.1))],'|','LineWidth',0.5,'Color','black');
text(((XC+XD)/2),(YB+(skala*0.12)),['\it L_{3} = ' num2str(L3*cosd(kat_prawa)*skala) ' cm'], 'Color', 'black', 'FontSize',11, 'HorizontalAlignment','center');
plot([XD XD], [YD (YB+(skala*0.1))], '-','LineWidth',0.1,'Color','black');
%LK
plot([0 XD],[(YB+(skala*0.1)) (YB+(skala*0.1))],'-','LineWidth',0.5,'Color','black');
plot([0 XD],[(YB+(skala*0.1)) (YB+(skala*0.1))],'|','LineWidth',0.5,'Color','black');
text((XD/2),(YB+(skala*0.12)),['\it L_{K} = ' num2str(LK*skala) ' cm'], 'Color', 'black', 'FontSize',11, 'HorizontalAlignment','center');
plot([0 0], [0 (YB+(skala*0.1))], '-','LineWidth',0.1,'Color','black');

%                   lewa strona

%LK
plot([0 XE],[(YB+(skala*0.1)) (YB+(skala*0.1))],'-','LineWidth',0.5,'Color','black');
plot([0 XE],[(YB+(skala*0.1)) (YB+(skala*0.1))],'|','LineWidth',0.5,'Color','black');
text((XE/2),(YB+(skala*0.12)),['\it 0.5 L_{K} = ' num2str(LK*100*0.5) ' cm'], 'Color', 'black', 'FontSize',11, 'HorizontalAlignment','center');
plot([XE XE], [YE (YB+(skala*0.1))], '-','LineWidth',0.1,'Color','black');
%L3
plot([XE XF],[(YB+(skala*0.1)) (YB+(skala*0.1))],'-','LineWidth',0.5,'Color','black');
plot([XE XF],[(YB+(skala*0.1)) (YB+(skala*0.1))],'|','LineWidth',0.5,'Color','black');
text(((XE+XF)/2),(YB+(skala*0.14)),['\it' num2str(cosd(kat_rozw)) ' L_{3} = ' num2str(L3*100*cosd(kat_rozw)) ' cm'], 'Color', 'black', 'FontSize',11, 'HorizontalAlignment','center');
plot([XF XF], [YF (YB+(skala*0.1))], '-','LineWidth',0.1,'Color','black');
%L2
plot([XF XG],[(YB+(skala*0.1)) (YB+(skala*0.1))],'-','LineWidth',0.5,'Color','black');
plot([XF XG],[(YB+(skala*0.1)) (YB+(skala*0.1))],'|','LineWidth',0.5,'Color','black');
text(((XF+XG)/2),(YB+(skala*0.16)),['\it' num2str(cosd(kat_rozw)) ' L_{2}cos(' num2str(theta_2) '\circ) = ' num2str(L2*cosd(theta_2)*100*cosd(kat_rozw)) ' cm'], 'Color', 'black', 'FontSize',11, 'HorizontalAlignment','center');
plot([XG XG], [YF (YB+(skala*0.1))], '-','LineWidth',0.1,'Color','black');
%L1
plot([XG XH],[(YB+(skala*0.1)) (YB+(skala*0.1))],'-','LineWidth',0.5,'Color','black');
plot([XG XH],[(YB+(skala*0.1)) (YB+(skala*0.1))],'|','LineWidth',0.5,'Color','black');
text(((XG+XH)/2),(YB+(skala*0.12)),['\it' num2str(cosd(kat_rozw)) ' L_{1}cos(' num2str(theta_1) '\circ) = ' num2str(L1*cosd(theta_1)*100*cosd(kat_rozw)) ' cm'], 'Color', 'black', 'FontSize',11, 'HorizontalAlignment','right');
plot([XH XH], [YH (YB+(skala*0.1))], '-','LineWidth',0.1,'Color','black');



% -------------- SIŁY -----------------



%wektor Q korpus
plot([0 0],[0 -(Q_korpus)], 'Color', 'c', 'LineWidth', 1.5 );
plot(0,-(Q_korpus), 'v','MarkerEdgeColor','c', 'MarkerFaceColor','c','MarkerSize',12);
text(0 + 0.005*skala,-(Q_korpus),['Q_{korpus} = ' num2str(Q_korpus) ' N'],'Color','cyan','FontSize',14, 'HorizontalAlignment','left');

%wektor Q_1
plot([XB XB],[YB (YB-(Q_1*5))], 'Color', 'c', 'LineWidth', 1.5 );
plot(XB,(YB-(Q_1*5)), 'v','MarkerEdgeColor','c', 'MarkerFaceColor','c','MarkerSize',12);
text(XB + 0.005*skala,(YB-(Q_1*5)),['Q_{1} = ' num2str(Q_1) ' N'],'Color','cyan','FontSize',14, 'HorizontalAlignment','left');

plot([XG XG],[YG (YG-(Q_1*5))*2], 'Color', 'c', 'LineWidth', 1.5 );
plot(XG,(YG-(Q_1*5))*2, 'v','MarkerEdgeColor','c', 'MarkerFaceColor','c','MarkerSize',12);
text(XG + 0.005*skala,(YG-(Q_1*5))*2,['2 Q_{1} = ' num2str(Q_1*2) ' N'],'Color','cyan','FontSize',14, 'HorizontalAlignment','left');

%wektor Q_2
plot([XC XC],[YC (YC-(Q_2*5))], 'Color', 'c', 'LineWidth', 1.5 );
plot(XC,(YC-(Q_2*5)), 'v','MarkerEdgeColor','c', 'MarkerFaceColor','c','MarkerSize',12);
text(XC + 0.005*skala,(YC-(Q_2*5)),['Q_{2} = ' num2str(Q_2) ' N'],'Color','cyan','FontSize',14, 'HorizontalAlignment','left');

plot([XF XF],[YF (YF-(Q_2*5))*2], 'Color', 'c', 'LineWidth', 1.5 );
plot(XF,(YF-(Q_2*5))*2, 'v','MarkerEdgeColor','c', 'MarkerFaceColor','c','MarkerSize',12);
text(XF - 0.005*skala,(YF-(Q_2*5))*2,['2 Q_{2} = ' num2str(Q_2*2) ' N'],'Color','cyan','FontSize',14, 'HorizontalAlignment','right');

%wektor Q_3
plot([XD XD],[YD (YD-(Q_3*5))], 'Color', 'c', 'LineWidth', 1.5 );
plot(XD,(YD-(Q_3*5)), 'v','MarkerEdgeColor','c', 'MarkerFaceColor','c','MarkerSize',12);
text(XD - 0.005*skala,(YD-(Q_3*5)),['Q_{3} = ' num2str(Q_3) ' N'],'Color','cyan','FontSize',14, 'HorizontalAlignment','right');

plot([XE XE],[YE (YE-(Q_3*5))*2], 'Color', 'c', 'LineWidth', 1.5 );
plot(XE,(YE-(Q_3*5))*2, 'v','MarkerEdgeColor','c', 'MarkerFaceColor','c','MarkerSize',12);
text(XE + 0.005*skala,(YE-(Q_3*5))*2,['2 Q_{3} = ' num2str(Q_3*2) ' N'],'Color','cyan','FontSize',14, 'HorizontalAlignment','left');

%wektor R1
plot([XA ((LK + (L3 + L2*cosd(theta_2) + 0.7*L1*cosd(theta_1))*cosd(kat_prawa) ) *skala)],[YA ((L2*sind(theta_2) - 0.7*L1*sind(theta_1)) *skala)], 'Color', 'r', 'LineWidth', 1.5 );
plot(((LK + (L3 + L2*cosd(theta_2) + 0.7*L1*cosd(theta_1))*cosd(kat_prawa) ) *skala),((L2*sind(theta_2) - 0.7*L1*sind(theta_1)) *skala),'^','MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',12);
text(((LK + (L3 + L2*cosd(theta_2) + 0.7*L1*cosd(theta_1))*cosd(kat_prawa) ) *skala) +0.01*skala,((L2*sind(theta_2) - 0.7*L1*sind(theta_1)) *skala),['R_{1} = ' num2str(R1) ' N'],'Color','red','FontSize',14, 'HorizontalAlignment','left');

%wektor R2
plot([XH -((LK*0.5 + L3*cosd(kat_rozw) + L2*cosd(theta_2)*cosd(kat_rozw) + 0.3*L1*cosd(theta_1)*cosd(kat_rozw)) *skala)],[YH ((L2*sind(theta_2) - 0.3*L1*sind(theta_1)) *skala)], 'Color', 'r', 'LineWidth', 1.5 );
plot(-((LK*0.5 + L3*cosd(kat_rozw) + L2*cosd(theta_2)*cosd(kat_rozw) + 0.3*L1*cosd(theta_1)*cosd(kat_rozw)) *skala),((L2*sind(theta_2) - 0.3*L1*sind(theta_1)) *skala),'^','MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',12);
text(-((LK*0.5 + L3*cosd(kat_rozw) + L2*cosd(theta_2)*cosd(kat_rozw) + 0.3*L1*cosd(theta_1)*cosd(kat_rozw)) *skala) -0.01*skala,((L2*sind(theta_2) - 0.3*L1*sind(theta_1)) *skala),['2 R_{2} = ' num2str(R2*2) ' N'],'Color','red','FontSize',14, 'HorizontalAlignment','right');

hold off;


% ------------------------ Rzut z gory --------------------------


figure();
grid on;
title(' Rzut z gory na robota. ', 'FontSize', 15);

line([0 LK*0.5*skala (LK*0.5*skala + LK*skala*cosd(60)) (LK*0.5*skala) -(LK*0.5*skala) -(LK*skala*cosd(60) + LK*0.5*skala) -LK*0.5*skala 0],[0 0 LK*skala*sind(60) 2*LK*skala*sind(60) 2*LK*skala*sind(60) LK*skala*sind(60) 0 0], 'Color', 'b', 'LineWidth', 1.5);
hold on;

%linia od prawej nogi - srodek

Y_A2 = LK*skala*sind(60);

p1 = line([(LK*skala) (LK*skala + L3*skala*cosd(kat_prawa)) (LK*skala + L3*skala*cosd(kat_prawa)+ L2*cosd(theta_2)*cosd(kat_prawa)*skala) (LK*skala + L3*skala*cosd(kat_prawa) + L2*skala*cosd(theta_2)*cosd(kat_prawa) + L1*skala*cosd(theta_1)*cosd(kat_prawa)) ],[Y_A2 (Y_A2 + L3*sind(kat_prawa)*skala) (Y_A2 + L2*sind(kat_prawa)*cosd(theta_2)*skala + L3*sind(kat_prawa)*skala) (Y_A2 +L3*sind(kat_prawa)*skala +L2*sind(kat_prawa)*cosd(theta_2)*skala +L1*sind(kat_prawa)*cosd(theta_1)*skala)], 'Color', 'b', 'LineWidth', 2);

plot((LK*skala + L3*skala*cosd(kat_prawa) + L2*skala*cosd(theta_2)*cosd(kat_prawa) + L1*skala*cosd(theta_1)*cosd(kat_prawa)),(Y_A2 +L3*sind(kat_prawa)*skala +L2*sind(kat_prawa)*cosd(theta_2)*skala +L1*sind(kat_prawa)*cosd(theta_1)*skala),'o', 'MarkerSize',10, 'MarkerEdgeColor','r', 'MarkerFaceColor','white');
text((LK*skala + L3*skala*cosd(kat_prawa) + L2*skala*cosd(theta_2)*cosd(kat_prawa) + L1*skala*cosd(theta_1)*cosd(kat_prawa)) -(0.001*skala),(Y_A2 +L3*sind(kat_prawa)*skala +L2*sind(kat_prawa)*cosd(theta_2)*skala +L1*sind(kat_prawa)*cosd(theta_1)*skala) -(0.005*skala),'A','Color','red','FontSize',12,  'HorizontalAlignment','left');

plot((LK*skala + L3*skala*cosd(kat_prawa) + L2*skala*cosd(theta_2)*cosd(kat_prawa)),(Y_A2 + L2*sind(kat_prawa)*cosd(theta_2)*skala + L3*sind(kat_prawa)*skala),'o', 'MarkerSize',10, 'MarkerEdgeColor','cyan', 'MarkerFaceColor','white');
text((LK*skala + L3*skala*cosd(kat_prawa) + L2*skala*cosd(theta_2)*cosd(kat_prawa)) -(0.001*skala),(Y_A2 + L2*sind(kat_prawa)*cosd(theta_2)*skala + L3*sind(kat_prawa)*skala) -(0.005*skala),'B','Color','cyan','FontSize',12,  'HorizontalAlignment','left');

plot((LK*skala + L3*skala*cosd(kat_prawa)),(Y_A2 + L3*sind(kat_prawa)*skala),'o', 'MarkerSize',10, 'MarkerEdgeColor','green', 'MarkerFaceColor','white');
text((LK*skala + L3*skala*cosd(kat_prawa)) -(0.001*skala),(Y_A2 + L3*sind(kat_prawa)*skala) -(0.005*skala),'C','Color','green','FontSize',12,  'HorizontalAlignment','left');

plot( (LK*skala),(Y_A2),'o', 'MarkerSize',10, 'MarkerEdgeColor','blue', 'MarkerFaceColor','white');
text((LK*skala) -(0.001*skala),(Y_A2) -(0.005*skala),'D','Color','blue','FontSize',12,  'HorizontalAlignment','left');

%linia od lewej nogi - srodek
line([-(LK*skala) -(LK*skala + L3*skala*cosd(kat_prawa)) -(LK*skala + L3*skala*cosd(kat_prawa)+ L2*cosd(theta_2)*cosd(kat_prawa)*skala) -(LK*skala + L3*skala*cosd(kat_prawa) + L2*skala*cosd(theta_2)*cosd(kat_prawa) + L1*skala*cosd(theta_1)*cosd(kat_prawa)) ],[Y_A2 (Y_A2 + L3*sind(kat_prawa)*skala) (Y_A2 + L2*sind(kat_prawa)*cosd(theta_2)*skala + L3*sind(kat_prawa)*skala) (Y_A2 +L3*sind(kat_prawa)*skala +L2*sind(kat_prawa)*cosd(theta_2)*skala +L1*sind(kat_prawa)*cosd(theta_1)*skala)], 'Color', 'g', 'LineWidth', 2, 'LineStyle','--');

plot(-(LK*skala + L3*skala*cosd(kat_prawa) + L2*skala*cosd(theta_2)*cosd(kat_prawa) + L1*skala*cosd(theta_1)*cosd(kat_prawa)),(Y_A2 +L3*sind(kat_prawa)*skala +L2*sind(kat_prawa)*cosd(theta_2)*skala +L1*sind(kat_prawa)*cosd(theta_1)*skala),'o', 'MarkerSize',10, 'MarkerEdgeColor','r', 'MarkerFaceColor','white');
text(-(LK*skala + L3*skala*cosd(kat_prawa) + L2*skala*cosd(theta_2)*cosd(kat_prawa) + L1*skala*cosd(theta_1)*cosd(kat_prawa)) -(0.001*skala),(Y_A2 +L3*sind(kat_prawa)*skala +L2*sind(kat_prawa)*cosd(theta_2)*skala +L1*sind(kat_prawa)*cosd(theta_1)*skala) -(0.005*skala),'H1','Color','red','FontSize',12,  'HorizontalAlignment','right');

plot(-(LK*skala + L3*skala*cosd(kat_prawa) + L2*skala*cosd(theta_2)*cosd(kat_prawa)),(Y_A2 + L2*sind(kat_prawa)*cosd(theta_2)*skala + L3*sind(kat_prawa)*skala),'o', 'MarkerSize',10, 'MarkerEdgeColor','cyan', 'MarkerFaceColor','white');
text(-(LK*skala + L3*skala*cosd(kat_prawa) + L2*skala*cosd(theta_2)*cosd(kat_prawa)) -(0.001*skala),(Y_A2 + L2*sind(kat_prawa)*cosd(theta_2)*skala + L3*sind(kat_prawa)*skala) -(0.005*skala),'G1','Color','cyan','FontSize',12,  'HorizontalAlignment','right');

plot(-(LK*skala + L3*skala*cosd(kat_prawa)),(Y_A2 + L3*sind(kat_prawa)*skala),'o', 'MarkerSize',10, 'MarkerEdgeColor','green', 'MarkerFaceColor','white');
text(-(LK*skala + L3*skala*cosd(kat_prawa)) -(0.001*skala),(Y_A2 + L3*sind(kat_prawa)*skala) -(0.005*skala),'F1','Color','green','FontSize',12,  'HorizontalAlignment','right');

plot(-(LK*skala),(Y_A2),'o', 'MarkerSize',10, 'MarkerEdgeColor','blue', 'MarkerFaceColor','white');
text(-(LK*skala) -(0.001*skala),(Y_A2) -(0.005*skala),'E1','Color','blue','FontSize',12,  'HorizontalAlignment','right');

% linia prawa noga - gora
Y_PG = 2*skala*LK*sind(60);
p2=line([LK*0.5*skala (LK*0.5*skala + L3*cosd(kat_rozw)*skala) ( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala) ( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala + L1*cosd(theta_1)*cosd(kat_rozw)*skala) ],[ Y_PG (Y_PG + L3*sind(kat_rozw)*skala) (Y_PG + L3*sind(kat_rozw)*skala +L2*cosd(theta_2)*sind(kat_rozw)*skala) (Y_PG + L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala + L1*cosd(theta_1)*sind(kat_rozw)*skala)],'Color', 'g', 'LineWidth', 2, "LineStyle",'--');

plot(( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala + L1*cosd(theta_1)*cosd(kat_rozw)*skala), (Y_PG + L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala + L1*cosd(theta_1)*sind(kat_rozw)*skala), 'o', 'MarkerSize',10, 'MarkerEdgeColor','r', 'MarkerFaceColor','white');
text(( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala + L1*cosd(theta_1)*cosd(kat_rozw)*skala),(Y_PG + L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala + L1*cosd(theta_1)*sind(kat_rozw)*skala) -(0.004*skala),'A1','Color','red','FontSize', 12,  'HorizontalAlignment','left');

plot(( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala),(Y_PG + L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala), 'o', 'MarkerSize',10, 'MarkerEdgeColor','cyan', 'MarkerFaceColor','white');
text(( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala),(Y_PG + L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala) -(0.004*skala),'B1','Color','cyan','FontSize', 12,  'HorizontalAlignment','left');

plot(( LK*0.5*skala + L3*cosd(kat_rozw)*skala),(Y_PG + L3*sind(kat_rozw)*skala), 'o', 'MarkerSize',10, 'MarkerEdgeColor','green','MarkerFaceColor','white');
text(( LK*0.5*skala + L3*cosd(kat_rozw)*skala),(Y_PG + L3*sind(kat_rozw)*skala) -(0.004*skala),'C1','Color','green','FontSize', 12,'HorizontalAlignment','left');

plot(( LK*0.5*skala),(Y_PG), 'o', 'MarkerSize',10, 'MarkerEdgeColor','blue','MarkerFaceColor','white');
text(( LK*0.5*skala),(Y_PG) -(0.004*skala),'D1','Color','blue','FontSize', 12,'HorizontalAlignment','left');

% linia lewa noga - gora

line([-LK*0.5*skala -(LK*0.5*skala + L3*cosd(kat_rozw)*skala) -( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala) -( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala + L1*cosd(theta_1)*cosd(kat_rozw)*skala) ],[ Y_PG (Y_PG + L3*sind(kat_rozw)*skala) (Y_PG + L3*sind(kat_rozw)*skala +L2*cosd(theta_2)*sind(kat_rozw)*skala) (Y_PG + L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala + L1*cosd(theta_1)*sind(kat_rozw)*skala)],'Color', 'b', 'LineWidth', 2, "LineStyle",'-');

plot(-( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala + L1*cosd(theta_1)*cosd(kat_rozw)*skala), (Y_PG + L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala + L1*cosd(theta_1)*sind(kat_rozw)*skala), 'o', 'MarkerSize',10, 'MarkerEdgeColor','r', 'MarkerFaceColor','white');
text(-( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala + L1*cosd(theta_1)*cosd(kat_rozw)*skala),(Y_PG + L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala + L1*cosd(theta_1)*sind(kat_rozw)*skala) -(0.004*skala),'H2','Color','red','FontSize', 12,  'HorizontalAlignment','right');

plot(-( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala),(Y_PG + L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala), 'o', 'MarkerSize',10, 'MarkerEdgeColor','cyan', 'MarkerFaceColor','white');
text(-( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala),(Y_PG + L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala) -(0.004*skala),'G2','Color','cyan','FontSize', 12,  'HorizontalAlignment','right');

plot(-( LK*0.5*skala + L3*cosd(kat_rozw)*skala),(Y_PG + L3*sind(kat_rozw)*skala), 'o', 'MarkerSize',10, 'MarkerEdgeColor','green','MarkerFaceColor','white');
text(-( LK*0.5*skala + L3*cosd(kat_rozw)*skala),(Y_PG + L3*sind(kat_rozw)*skala) -(0.004*skala),'F2','Color','green','FontSize', 12,'HorizontalAlignment','right');

plot(-( LK*0.5*skala),(Y_PG), 'o', 'MarkerSize',10, 'MarkerEdgeColor','blue','MarkerFaceColor','white');
text(-( LK*0.5*skala),(Y_PG) -(0.004*skala),'E2','Color','blue','FontSize', 12,'HorizontalAlignment','right');

% linia prawa noga - dol
line([ LK*0.5*skala (LK*0.5*skala + L3*cosd(kat_rozw)*skala) ( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala) ( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala + L1*cosd(theta_1)*cosd(kat_rozw)*skala) ],[ 0 (-L3*sind(kat_rozw)*skala) -(L3*sind(kat_rozw)*skala +L2*cosd(theta_2)*sind(kat_rozw)*skala) -(L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala + L1*cosd(theta_1)*sind(kat_rozw)*skala)],'Color', 'g', 'LineWidth', 2, "LineStyle",'--');

plot(( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala + L1*cosd(theta_1)*cosd(kat_rozw)*skala), -(L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala + L1*cosd(theta_1)*sind(kat_rozw)*skala), 'o', 'MarkerSize',10, 'MarkerEdgeColor','r', 'MarkerFaceColor','white');
text(( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala + L1*cosd(theta_1)*cosd(kat_rozw)*skala), -(L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala + L1*cosd(theta_1)*sind(kat_rozw)*skala) -(0.004*skala),'A2','Color','red','FontSize', 12,  'HorizontalAlignment','right');

plot(( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala),-(L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala), 'o', 'MarkerSize',10, 'MarkerEdgeColor','cyan', 'MarkerFaceColor','white');
text(( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala),-(L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala) -(0.004*skala),'B2','Color','cyan','FontSize', 12,  'HorizontalAlignment','right');

plot(( LK*0.5*skala + L3*cosd(kat_rozw)*skala),-(L3*sind(kat_rozw)*skala), 'o', 'MarkerSize',10, 'MarkerEdgeColor','green','MarkerFaceColor','white');
text(( LK*0.5*skala + L3*cosd(kat_rozw)*skala),-(L3*sind(kat_rozw)*skala) -(0.004*skala),'C2','Color','green','FontSize', 12,'HorizontalAlignment','right');

plot(( LK*0.5*skala),0, 'o', 'MarkerSize',10, 'MarkerEdgeColor','blue','MarkerFaceColor','white');
text(( LK*0.5*skala), -(0.004*skala),'D2','Color','blue','FontSize', 12,'HorizontalAlignment','right');

% linia lewa noga - dol

if( kroczenie == 1)
    kat_krocz = -kat_rozw;
else
    kat_krocz = kat_rozw;
end    

line([ -LK*0.5*skala -(LK*0.5*skala + L3*cosd(kat_krocz)*skala) -( LK*0.5*skala + L3*cosd(kat_krocz)*skala + L2*cosd(theta_2)*cosd(kat_krocz)*skala) -( LK*0.5*skala + L3*cosd(kat_krocz)*skala + L2*cosd(theta_2)*cosd(kat_krocz)*skala + L1*cosd(theta_1)*cosd(kat_krocz)*skala) ],[ 0 (-L3*sind(kat_krocz)*skala) -(L3*sind(kat_krocz)*skala +L2*cosd(theta_2)*sind(kat_krocz)*skala) -(L3*sind(kat_krocz)*skala + L2*cosd(theta_2)*sind(kat_krocz)*skala + L1*cosd(theta_1)*sind(kat_krocz)*skala)],'Color', 'b', 'LineWidth', 2, "LineStyle",'-');

plot(-( LK*0.5*skala + L3*cosd(kat_krocz)*skala + L2*cosd(theta_2)*cosd(kat_krocz)*skala + L1*cosd(theta_1)*cosd(kat_krocz)*skala), -(L3*sind(kat_krocz)*skala + L2*cosd(theta_2)*sind(kat_krocz)*skala + L1*cosd(theta_1)*sind(kat_krocz)*skala), 'o', 'MarkerSize',10, 'MarkerEdgeColor','r', 'MarkerFaceColor','white');
text(-( LK*0.5*skala + L3*cosd(kat_krocz)*skala + L2*cosd(theta_2)*cosd(kat_krocz)*skala + L1*cosd(theta_1)*cosd(kat_krocz)*skala), -(L3*sind(kat_krocz)*skala + L2*cosd(theta_2)*sind(kat_krocz)*skala + L1*cosd(theta_1)*sind(kat_krocz)*skala) -(0.004*skala),'H','Color','red','FontSize', 12,  'HorizontalAlignment','left');

plot(-( LK*0.5*skala + L3*cosd(kat_krocz)*skala + L2*cosd(theta_2)*cosd(kat_krocz)*skala),-(L3*sind(kat_krocz)*skala + L2*cosd(theta_2)*sind(kat_krocz)*skala), 'o', 'MarkerSize',10, 'MarkerEdgeColor','cyan', 'MarkerFaceColor','white');
text(-( LK*0.5*skala + L3*cosd(kat_krocz)*skala + L2*cosd(theta_2)*cosd(kat_krocz)*skala),-(L3*sind(kat_krocz)*skala + L2*cosd(theta_2)*sind(kat_krocz)*skala) -(0.004*skala),'G','Color','cyan','FontSize', 12,  'HorizontalAlignment','left');

plot(-( LK*0.5*skala + L3*cosd(kat_krocz)*skala),-(L3*sind(kat_krocz)*skala), 'o', 'MarkerSize',10, 'MarkerEdgeColor','green','MarkerFaceColor','white');
text(-( LK*0.5*skala + L3*cosd(kat_krocz)*skala),-(L3*sind(kat_krocz)*skala) -(0.004*skala),'F','Color','green','FontSize', 12,'HorizontalAlignment','left');

plot(-( LK*0.5*skala),0, 'o', 'MarkerSize',10, 'MarkerEdgeColor','blue','MarkerFaceColor','white');
text(-( LK*0.5*skala), -(0.004*skala),'E','Color','blue','FontSize', 12,'HorizontalAlignment','left');

% srodek ciezkosci
plot(0, Y_A2, 'square', 'MarkerSize', 15, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white');
text(0, (Y_A2 - 0.01*skala), ' śr. ciężkości korpusu ','Color','black','FontSize', 12,'HorizontalAlignment', 'center' );



% ----------- katy -------------

%definiuj kat rozwarcia - 1ewa strona
rozdz_thrz = ( 180 - kat_rozw_pocz):0.1:180;

%centrum kata
r_rz = skala*0.03;

% def rownania
X_rz = -( LK * 0.5 * skala ) + r_rz*cosd(rozdz_thrz);
Y_rz = Y_PG + r_rz*sind(rozdz_thrz);

plot(X_rz,Y_rz, 'Color', 'black','LineStyle','-', 'LineWidth', 0.5);
if( kat_rozw_pocz > 10 )
    text((-( LK * 0.5 * skala )*2.25),(Y_PG +(skala*0.005)),['\it \theta_{rozwarcie} = ' num2str(kat_rozw) '\circ'], 'Color', 'black', 'FontSize',13, 'HorizontalAlignment','left');
elseif ( kat_rozw_pocz >= 0 && kat_rozw_pocz <= 10)
    text((-( LK * 0.5 * skala )*2.25),(Y_PG -(skala*0.01)),['\it \theta_{rozwarcie} = ' num2str(kat_rozw) '\circ'], 'Color', 'black', 'FontSize',13, 'HorizontalAlignment','left');
else 
    text((-( LK * 0.5 * skala )*2.25),(Y_PG -(skala*0.01)),['\it \theta_{rozwarcie} = ' num2str(kat_rozw) '\circ'], 'Color', 'black', 'FontSize',13, 'HorizontalAlignment','left');
end    

line([ -(LK * 0.5 * skala) (X_rz - skala*0.01)],[Y_PG Y_PG], 'Color', 'Black', 'LineStyle', '--', 'LineWidth', 0.5);


% ---------- kat prawa ------------------
%definiuj kat rozwarcia dla srodkowej nogi - prawa strona
prawa_thrz = 0:0.1:kat_prawa;

% def rownania prawej
X_pz = LK * skala + r_rz*cosd(prawa_thrz);
Y_pz = Y_A2 + r_rz*sind(prawa_thrz);

plot(X_pz,Y_pz, 'Color', 'black','LineStyle','-', 'LineWidth', 0.5);
if( kat_prawa > 10 )
    text((( LK * skala )*1.5),(Y_A2 + (skala*0.005)),['\it \theta_{prawa} = ' num2str(kat_prawa) '\circ'], 'Color', 'black', 'FontSize',13, 'HorizontalAlignment','right');
elseif ( kat_prawa >= 0 && kat_prawa <= 10)
    text((( LK * skala )*1.5),(Y_A2 - (skala*0.01)),['\it \theta_{prawa} = ' num2str(kat_prawa) '\circ'], 'Color', 'black', 'FontSize',13, 'HorizontalAlignment','right');
else 
    text((( LK * skala )*1.5),(Y_A2 - (skala*0.01)),['\it \theta_{prawa} = ' num2str(kat_prawa) '\circ'], 'Color', 'black', 'FontSize',13, 'HorizontalAlignment','right');
end    

line([ (LK * skala) (X_pz(end) + skala*0.01)],[Y_A2  Y_A2 ], 'Color', 'Black', 'LineStyle', '--', 'LineWidth', 0.5);

% granice osi y
ylim([((-(Y_PG + L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala + L1*cosd(theta_1)*sind(kat_rozw)*skala) + 14)- skala*0.025) ((Y_PG + L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala + L1*cosd(theta_1)*sind(kat_rozw)*skala) + skala*0.025)]);


% wielobok podparcia oceniajacy stabilnosc statyczna.
p3 = line([-( LK*0.5*skala + L3*cosd(kat_krocz)*skala + L2*cosd(theta_2)*cosd(kat_krocz)*skala + L1*cosd(theta_1)*cosd(kat_krocz)*skala) (LK*skala + L3*skala*cosd(kat_prawa) + L2*skala*cosd(theta_2)*cosd(kat_prawa) + L1*skala*cosd(theta_1)*cosd(kat_prawa)) -( LK*0.5*skala + L3*cosd(kat_rozw)*skala + L2*cosd(theta_2)*cosd(kat_rozw)*skala + L1*cosd(theta_1)*cosd(kat_rozw)*skala) -( LK*0.5*skala + L3*cosd(kat_krocz)*skala + L2*cosd(theta_2)*cosd(kat_krocz)*skala + L1*cosd(theta_1)*cosd(kat_krocz)*skala)],[-(L3*sind(kat_krocz)*skala + L2*cosd(theta_2)*sind(kat_krocz)*skala + L1*cosd(theta_1)*sind(kat_krocz)*skala) (Y_A2 +L3*sind(kat_prawa)*skala +L2*sind(kat_prawa)*cosd(theta_2)*skala +L1*sind(kat_prawa)*cosd(theta_1)*skala) (Y_PG + L3*sind(kat_rozw)*skala + L2*cosd(theta_2)*sind(kat_rozw)*skala + L1*cosd(theta_1)*sind(kat_rozw)*skala) -(L3*sind(kat_krocz)*skala + L2*cosd(theta_2)*sind(kat_krocz)*skala + L1*cosd(theta_1)*sind(kat_krocz)*skala)], 'Color', 'red', 'LineWidth', 1);


hold off;
legend([p1 p2 p3],{' Noga podpierajaca ',' Noga w powietrzu ', ' wielobok podparcia '});

end