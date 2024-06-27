%% Plot full EOM simulation and Slow Invariant manifold
% This code was used to generate the time plots in https://doi.org/10.1016/j.ymssp.2024.111615
% From Github repo: https://github.com/kevindekemele/PZT_cubic_shunt_NES 
% Author: Kevin Dekemele
%% Define system
%close all
clear all

w_i = 1;

zeta_e = 0.20;
r_i = 0.20;
k_i = 0.20;
zeta_i = 0.0073;
Fex1 = 0.41;

%% Time simulation

t = 0:0.1:3000;

%Initial conditions, should be non zero to reach isoalted curves
x0 = [0 0];
x0_dot = [0 0];

%forcing frequency
Om_f = 1;

% Equations of motion Eq. 8
f1 = @(t,y)[y(3);y(4);...
    -(y(1) + 2*zeta_i*y(3) + k_i*r_i*y(2) - Fex1(1)*cos(Om_f*t));...
    -( r_i^2*y(2) + 2*zeta_e*r_i*y(4) + k_i*r_i*y(1)+r_i*(k_i*y(1)+r_i*y(2))^3);
    ];

Prec = 1e-8;
    % Actual numerical simulation
    
options = odeset('RelTol',Prec,'AbsTol',[Prec Prec Prec Prec]);
[T1,Y1] = ode45(f1,t,[x0 x0_dot],options);

% Eq 8c
V = (k_i*Y1(:,1) + r_i*Y1(:,2))*r_i;
V_spd = (k_i*Y1(:,3) + r_i*Y1(:,4))*r_i;

% Envelope extraction, this has to be a bit tuned for each simulation as QP
% may have other beating frequencies/depth
envY  = envelope(abs(Y1(:,1)),80,'peak');
envB  = envelope(abs(V),130,'peak');

%% Plotting
figure
box on
subplot(2,1,1);
plot(T1,Y1(:,1),'k','LineWidth',2)
    hold on
    plot(T1,envY,'LineWidth',2)
    % numerical_SIM
  xlim([1000, 1300])
   ax = gca; 
ax.FontSize = 15; 
ylabel('$a$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis

subplot(2,1,2);
box on
    box on
hold on
    plot(T1,V,'k','LineWidth',2)
    hold on
    plot(T1,envB,'LineWidth',2)
    %xlim([500, 2000])
  xlim([1000, 1300])
   ax = gca; 
ax.FontSize = 15; 
ylabel('$b$','FontSize',16.5,'interpreter','latex') 
xlabel('$\tau$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
     

%% Plotting of SIM
      Zb = 0:0.001:1.8^2;  
      Za=((2*zeta_e*r_i*Om_f)^2+ ( (Om_f^2-r_i^2)-3/4*Zb).^2).*Zb./(r_i^2*k_i^2*(Om_f^4+(2*zeta_e*r_i*Om_f)^2));

    Om = linspace(0.85,1.15,100);
% Stability of SIM
  for j=2:length(Zb)
              
        alphaminbeta_SIM(j) = -acos(-((r_i^2-Om_f^2)*sqrt(Zb(j))-4*zeta_e^2*r_i^2*sqrt(Zb(j))+3/4*sqrt(Zb(j))^3)/(Om_f^2*sqrt(Za(j))*k_i*r_i+4*zeta_e^2*r_i^3*k_i*sqrt(Za(j))));
          
        alpha1_SIM(j) =  0;%-asin((2*zeta_i*sqrt(Za(j))*Om_f+sin(-alphaminbeta_SIM(j))*k_i/r_i*sqrt(Zb(j)))/Fex1);%-acos((1-k_i^2-Om(i)^2)*sqrt(zo(i,1))+cos(-alphaminbeta(i,1))*k_i/r_i*sqrt(zna(i,1))/Fex1);
        beta1_SIM(j) =   -(alphaminbeta_SIM(j)  - alpha1_SIM(j) );

          M= [
             -(r_i^2-Om_f^2)-2*zeta_e*1i*r_i*Om_f-6/4*Zb(j),  -3/4*Zb(j)*exp(1i*alphaminbeta_SIM(j))^2;
             -conj( -3/4*Zb(j)*exp(1i*alphaminbeta_SIM(j))^2), -conj(-(r_i^2-Om_f^2)-2*zeta_e*1i*r_i*Om_f-6/4*Zb(j));]/(Om_f*2*1i);
          eigM(j,:) =   eig(M);
         if(any(real(eig(M))>0))
               B_SIM_unstable(j) =sqrt(Zb(j)); 
               A_SIM_unstable(j) =sqrt(Za(j));
               B_SIM_stable(j) =NaN; 
               A_SIM_stable(j) =NaN;
            else
               B_SIM_stable(j) =sqrt(Zb(j));  
               A_SIM_stable(j) =sqrt(Za(j));  
               B_SIM_unstable(j) = NaN; 
               A_SIM_unstable(j) = NaN;
            end
          
               
  end


     figure(501)
         box on

     hold on
         plot(B_SIM_unstable,A_SIM_unstable,'k--','LineWidth',2)
    plot(B_SIM_stable,A_SIM_stable,'k','LineWidth',2)
 ax = gca; 
ax.FontSize = 15; 
xlabel('$b$','FontSize',16.5,'interpreter','latex') 
ylabel('$a$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
    plot(abs(envB(round(3*length(Y1)/4):round(3.9*length(Y1)/4))),abs(envY(round(3*length(Y1)/4):round(3.9*length(Y1)/4))),'LineWidth',2)
    
    %%

