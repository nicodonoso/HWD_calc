close all; clear all; 

%% INTEGRACIÓN NUMERICA 
%Constantes
%Vp     = 1.26e-2;%[m3s-1] from Greenler(2014) % flow into the hole CECs lake is 210 L/min = 3.5 L/s
Vp     = 0.35e-2;% [m3s-1] Vp BAS is 210 L/min -> 3.5 L/s =  0.0035 m3/s = 0.35e-2 (revisado)
Ttip   = 93;%[C] temperature of water at point of drill system
%v     = 0.0375; %[m s-1] velocity of drill, CECs lake is 0.7 - 2.0 m/min --> 0.0117-0.0333 m/s
v      = [0.0117;0.0333];%[m s-1] velocity of drill in between 0.0333-0.0117 m/s 
%Tinf  = -50; %[C] Temperature of ice
Tinf   = -27; %[C] Temperature of ice CECs (surface)
rho_i  = 917; %[kg m-3] Density of ice
rho_w  = 982; %[kg m-3] Density of water
c_w    = 4170; % [J kg-1 c-1] specific heat of water, Greenler(2014)
c_i    = 1950; % [J kg-1 c-1] specific heat of ice, Greenler(2014)
c_f    = 335000;% [J kg-1] Heat of fusion of ice, Greenler(2014)
Dm     = 0.93;% Dimensionless, ratio of specific volume of water to that ice, Greenler(2014)
k_w    = 0.655;% [W m-1 C-1] % conductividad termica del agua?
k_i    = 2.2;% [W m-1C-1] Thermal conductivity of ice, Greenler(2014)
k_h    = 0.37;% [W m-1C-1] Thermal conductivity of house, enviado por Jano hwatsapp miércoles 29 abril
pi     = 3.1415;% Babylonian tablet (ca. 1900–1680 BC)
r1     = 0.0159;%[m] 15.9 mm internal radio of hose
r2     = 0.0221;%[m] 22.1 mm external radio of hose
R      = 0.125;%[m] 25 cm diameter of mean hole to drill 
lake_depth = 2650;%[m] from CECs team google spread sheet 

Z_h    = log(r2/r1)/(2*pi*k_h);% Efective thermnal resistence of the hose
Z_d    = log(R/r2)/(2*pi*k_w);% Efective thermnal resistence of the water
Labda  = rho_w*c_w*Vp*(Z_h+Z_d);% Characteristic decay length, assumed to be independent of depth. Humphrey and Echelmeyer (1990) ec(2) 
% Abreviaciones:
A = 2*v*rho_i*(c_f-c_i*Tinf)/(k_w*0.023);%revisar si es k_w el k
B = (rho_i*(c_f-c_i*Tinf))/(rho_w*c_w);
C = v(1)*rho_i*(c_f-c_i*Tinf);


% %% Evaluamos ecuación (1) de Greenler(2014)
R = 0.0001:0.001:0.7;
for j=1:length(v)
    for i=1:length(R)
        T_w(i,j) = (Ttip*Vp)/(Vp+Dm*pi*v(j)*R(i)*R(i))-(pi*v(j)*R(i)*R(i))*B/(Vp+Dm*pi*v(j)*R(i)*R(i));
    end
end
ff1 = figure(1) %Greenler(2014)
plot(R,T_w(:,1),'LineWidth',4)
hold on
plot(R,T_w(:,2),'LineWidth',4)
xlabel('Hole radio [m]')
ylabel('Pump out water temperature [ºC]')
title(['T_{tip} = ',num2str(Ttip),'{\circ}C, Flow = ',num2str(Vp*60*1000),' L min^{-1}'])
legend(['v_{drill} = ',num2str(v(1)*60-0.002),' [m min^{-1}]'], ['v_{drill} = ',num2str(v(2)*60+0.002),' [m min^{-1}]'])
grid on
%axis([0.05 0.3 0 80])
ylim([0 Ttip*1.1])
%disp(['La temperatura del agua es de: ',num2str(T_w),' para un radio de ',num2str(R)])
%print(ff1,'-djpeg','Greenler_et_al_ec1','-r600')

%% Radio máximo variando flujo
% if T_w = 0 mean all heat is realese so we get the max radio
% 1 How Rmax change in change of flow into the hole
clear Vp
Vp = 0:0.001e-2:0.7e-2; 
for j=1:length(v)
    for i = 1:length(Vp)% Greenler(2014) Ec(3) 

        Rmax(i,j) = sqrt( (Vp(i)*Ttip*rho_w*c_w)/ (pi*v(j)*rho_i*(c_f-c_i*Tinf)) );
    
    end
end

ff2 = figure(2)
plot(Vp*1000,Rmax(:,1),'LineWidth',4)
hold on
plot(Vp*1000,Rmax(:,2),'LineWidth',4)
plot(3.5,0.125,'*r')
legend('Theory, max. efficiency v_{Drill}=0.7 m/min','Theory, max. efficiency v_{Drill}=2.0 m/min','BAS Drill(3.5 L/s,12.5 [cm])','Location','northwest')
xlabel('Water flow rate throug drill [L/s]')
ylabel('Max. radius [m]')
title(['T_{tip} = ',num2str(Ttip),' {\circ}C'])
grid on
%print(ff2,'-djpeg','Greenler_ec3_Drill_Flow','-r600')

%% Radio máximo variando velocidad de perforación
%Vp BAS is 210 L/min -> 3.5 L/s =  0.0035 m3/s = 0.35e-2 (revisado)

Vp = 0.35e-2; 
clear v Rmax 
v = 0.001:0.001:0.09;
for i=1:length(v)
   
     Rmax(i) = sqrt( (Vp*Ttip*rho_w*c_w)/ (pi*v(i)*rho_i*(c_f-c_i*Tinf)) );
    
end

ff3 = figure(3)
plot(v,Rmax,'LineWidth',4)
xlabel('Drill velocity [m s^{-1}]')
ylabel('Max. radius [m]')
title(['Flow = ',num2str(Vp*1000*60),'[L/min], T_{tip} = ',num2str(Ttip),'{\circ}C'])
grid on
%print(ff3,'-djpeg','Greenler_ec3_drill_velocity_2','-r600')
%% Calculo de la caida de Tº del flujo de perforación

y = 0:0.1:lake_depth;%[m]
for depth=1:length(y)
T_d(depth) = Ttip*exp(-y(depth)/Labda);% Humphrey and Echelmeyer (1990) ec(1) 
end

ff4 = figure(4)
plot(y,T_d,'LineWidth',4)
xlabel('Depth [m]')
ylabel('Nozzle tip temperature dropt [{\circ}C]')
grid on
xlim([0 lake_depth])
%print(ff4,'-djpeg','Temperature_drop_nozzle','-r600')

return
%% Función de forma del agujero

%FUN1 = T_w
fun1 = @(R)((Ttip*Vp)./(Vp+Dm*pi*v*R.^2)-(pi*v*R.^2)*B./(Vp+Dm*pi*v*R.^2));
% FUN2 = dy/dR
fun2 = @(R)(A*(0.00493*fun1(R)+0.055).^0.3/fun1(R)*(2*rho_w*v*R.*(27*fun1(R)+500).^0.8));

q = integral(fun2,0,0.3)%0.05)
