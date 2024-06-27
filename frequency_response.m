%% Plot dimensionless frequency response
% This code was used to generate the frequency responses in https://doi.org/10.1016/j.ymssp.2024.111615
% From Github repo: https://github.com/kevindekemele/PZT_cubic_shunt_NES 
% Author: Kevin Dekemele
%% Define system
clear all
%close all
w_i = 1;

zeta_e = 0.20;
r_i = 0.20;
k_i = 0.20;
zeta_i = 0.0073;
Fex1 =0.41;
%% FR computation
Om = linspace(0.85,1.15,5000);

for i=1:length(Om)
    
    %Eq.19
    a = ( (1-k_i^2-Om(i)^2)*3/4)^2+(2*zeta_i*Om(i)*3/4)^2;
    b = (1-k_i^2-Om(i)^2)*3/2*((1-k_i^2-Om(i)^2)*(r_i^2-Om(i)^2)-4*zeta_i*zeta_e*Om(i)^2*r_i-k_i^2*Om(i)^2) + zeta_i*Om(i)*3*( (1-Om(i)^2)*(Om(i)*2*zeta_e*r_i) + (2*zeta_i*Om(i))*(r_i^2-Om(i)^2) );
    c = ((1-k_i^2-Om(i)^2)*(r_i^2-Om(i)^2)-4*zeta_i*zeta_e*Om(i)^2*r_i-k_i^2*Om(i)^2)^2+ ( (1-Om(i)^2)*(Om(i)*2*zeta_e*r_i) + (2*zeta_i*Om(i))*(r_i^2-Om(i)^2) )^2;
    d = -Fex1(1)^2*( (2*Om(i)*zeta_e*r_i)^2+Om(i)^4)*(r_i*k_i)^2;
    [r] = roots([a,b,c,d]);
    
    % Open circuit FR
     zo_OC(i) = (-Om(i)^2 + 1i*Om(i)*2*zeta_i + (1))\Fex1;

    %Check if one or three solutions
    
    if(~isreal(r))
        r_real=r(imag(r)==0);
        %Eq. 17
        zo(i,1) = ((2*zeta_e*r_i*Om(i))^2+ ( (Om(i)^2-r_i^2)-3/4*r_real)^2)*r_real/(r_i^2*k_i^2*(Om(i)^4+(2*zeta_e*r_i*Om(i))^2));
        zo(i,2)=NaN;
        zo(i,3)=NaN;
        zna(i,1) = r_real;
        zna(i,2)=NaN;
        zna(i,3)=NaN;
        
        %Phase, computed from Eq. 20
        alphaminbeta(i,1) = acos(-((r_i^2-Om(i)^2)*sqrt(zna(i,1))-4*zeta_e^2*r_i^2*sqrt(zna(i,1))+3/4*sqrt(zna(i,1))^3)/(Om(i)^2*sqrt(zo(i,1))*k_i*r_i+4*zeta_e^2*r_i^3*k_i*sqrt(zo(i,1)))); 
        alpha1(i,1) =  -asin((2*zeta_i*sqrt(zo(i,1))*Om(i)+sin(-alphaminbeta(i,1))*k_i/r_i*sqrt(zna(i,1)))/Fex1);%-acos((1-k_i^2-Om(i)^2)*sqrt(zo(i,1))+cos(-alphaminbeta(i,1))*k_i/r_i*sqrt(zna(i,1))/Fex1);
        beta1(i,1) =   -(alphaminbeta(i,1)  - alpha1(i,1) );

        %Jacobian to compute stability Appendix A
        % The evolution of the eigenvalues could be saved for bifurcation
        % study
          M2= [-(1-k_i^2-Om(i)^2)-2*1i*zeta_i*Om(i), 0, -k_i/r_i, 0;
            0,  (1-k_i^2-Om(i)^2)-2*1i*zeta_i*Om(i),0,  k_i/r_i;
            2*zeta_e*1i*r_i^2*Om(i)*k_i-r_i*k_i*((1-k_i^2) + 2*1i*Om(i)*zeta_i), 0, -(r_i^2-Om(i)^2)-2*zeta_e*1i*r_i*Om(i)-6/4*zna(i,1)-r_i*k_i*k_i/r_i,  -3/4*zna(i,1)*exp(1i*beta1(i,1))^2;
            0,  -conj( 2*zeta_e*1i*r_i^2*Om(i)*k_i-r_i*k_i*((1-k_i^2) + 2*1i*Om(i)*zeta_i)),  -conj( -3/4*zna(i,1)*exp(1i*beta1(i,1))^2), -conj(-(r_i^2-Om(i)^2)-2*zeta_e*1i*r_i*Om(i)-6/4*zna(i,1)-r_i*k_i*k_i/r_i);]/(Om(i)*2*1i );
       EIGM2(i,:) =  sort(eig(M2));
       
       %determine which a and b are stable/unstable
            if(any(real(eig(M2))>0))
               B_unstable(i,1) =sqrt(zna(i,1)); 
               A_unstable(i,1) =sqrt(zo(i,1)); 
               B_stable(i,1) = NaN; 
               A_stable(i,1) = NaN;
            else
           B_stable(i,1) =sqrt(zna(i,1)); 
               A_stable(i,1) =sqrt(zo(i,1));  
               B_unstable(i,1) = NaN; 
               A_unstable(i,1) = NaN;
            end
                 B_stable(i,2) = NaN; 
               A_stable(i,2) = NaN;
                    B_stable(i,3) = NaN; 
               A_stable(i,3) = NaN;
            B_unstable(i,2) = NaN; 
               A_unstable(i,2) = NaN;
               B_unstable(i,3) = NaN; 
               A_unstable(i,3) = NaN;
    else
            %sorting needed to gain nice plots in the end, perhaps there is
            %a better way
            r=sort(r,'ascend') ;

         zo(i,1)=((2*zeta_e*r_i*Om(i))^2+ ( (Om(i)^2-r_i^2)-3/4*r(1))^2)*r(1)/(r_i^2*k_i^2*(Om(i)^4+(2*zeta_e*r_i*Om(i))^2));
        zo(i,2)=((2*zeta_e*r_i*Om(i))^2+ ( (Om(i)^2-r_i^2)-3/4*r(2))^2)*r(2)/(r_i^2*k_i^2*(Om(i)^4+(2*zeta_e*r_i*Om(i))^2));
        zo(i,3)=((2*zeta_e*r_i*Om(i))^2+ ( (Om(i)^2-r_i^2)-3/4*r(3))^2)*r(3)/(r_i^2*k_i^2*(Om(i)^4+(2*zeta_e*r_i*Om(i))^2));

        zna(i,1) = r(1);
        zna(i,2)= r(2);
        zna(i,3)=r(3);


                for j=1:3
                                     alphaminbeta(i,j) = acos(-((r_i^2-Om(i)^2)*sqrt(zna(i,j))-4*zeta_e^2*r_i^2*sqrt(zna(i,j))+3/4*sqrt(zna(i,j))^3)/(Om(i)^2*sqrt(zo(i,j))*k_i*r_i+4*zeta_e^2*r_i^3*k_i*sqrt(zo(i,j))));
          
        alpha1(i,j) =  -asin((2*zeta_i*sqrt(zo(i,j))*Om(i)+sin(-alphaminbeta(i,j))*k_i/r_i*sqrt(zna(i,j)))/Fex1);%-acos((1-k_i^2-Om(i)^2)*sqrt(zo(i,1))+cos(-alphaminbeta(i,1))*k_i/r_i*sqrt(zna(i,1))/Fex1);
        beta1(i,j) =   -(alphaminbeta(i,j) - alpha1(i,1) );

          M2= [-(1-k_i^2-Om(i)^2)-2*1i*zeta_i*Om(i), 0, -k_i/r_i, 0;
            0,  (1-k_i^2-Om(i)^2)-2*1i*zeta_i*Om(i),0,  k_i/r_i;
            2*zeta_e*1i*r_i^2*Om(i)*k_i-r_i*k_i*((1-k_i^2) + 2*1i*Om(i)*zeta_i), 0, -(r_i^2-Om(i)^2)-2*zeta_e*1i*r_i*Om(i)-6/4*zna(i,j)-r_i*k_i*k_i/r_i,  -3/4*zna(i,j)*exp(1i*beta1(i,j))^2;
            0,  -conj( 2*zeta_e*1i*r_i^2*Om(i)*k_i-r_i*k_i*((1-k_i^2) + 2*1i*Om(i)*zeta_i)),  -conj( -3/4*zna(i,j)*exp(1i*beta1(i,j))^2), -conj(-(r_i^2-Om(i)^2)-2*zeta_e*1i*r_i*Om(i)-6/4*zna(i,j)-r_i*k_i*k_i/r_i);]/(Om(i)*2*1i );
    
        eig(M2)
            if(any(real(eig(M2))>0))
               B_unstable(i,j) =sqrt(zna(i,j)); 
               A_unstable(i,j) =sqrt(zo(i,j)); 
               B_stable(i,j) = NaN; 
               A_stable(i,j) = NaN;
            else
           B_stable(i,j) =sqrt(zna(i,j)); 
               A_stable(i,j) =sqrt(zo(i,j));  
               B_unstable(i,j) = NaN; 
               A_unstable(i,j) = NaN;
            end
                end
        
    end
    
    % maximum (=saturation amplitude) and minmum
    a2 = 27;
    b2= (48*r_i^2-48*Om(i)^2);
    c2 = 64*Om(i)^2*r_i^2*zeta_e^2+16*r_i^4-32*Om(i)^2*r_i^2+16*Om(i)^4;
    [r] = roots([a2,b2,c2]);
 if(isreal(r))
         r_min=min(r);
         zna_min(i) = r_min;
          zo_max(i) = ((2*zeta_e*r_i*Om(i))^2+ ( (Om(i)^2-r_i^2)-3/4*r_min)^2)*r_min/(r_i^2*k_i^2*(Om(i)^4+(2*zeta_e*r_i*Om(i))^2));
         r_max=max(r);
          zna_max(i) = r_max;
          zo_min(i) = ((2*zeta_e*r_i*Om(i))^2+ ( (Om(i)^2-r_i^2)-3/4*r_max)^2)*r_max/(r_i^2*k_i^2*(Om(i)^4+(2*zeta_e*r_i*Om(i))^2));
 end
            
end
%% Plotting
figure(1)
hold on
hold on
plot(Om,A_stable,'k','LineWidth',2)
hold on
plot(Om,A_unstable,'k--','LineWidth',2)
plot(Om,sqrt(zo_max),'k--')
plot(Om,sqrt(zo_min),'k--')
plot(Om,abs(zo_OC),'k:')
 ax = gca;
ax.FontSize = 15; 
xlabel('$\Omega$','FontSize',16.5,'interpreter','latex') 
ylabel('$a$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
figure(2)
hold on
plot(Om,B_stable,'k','LineWidth',2)
hold on
plot(Om,B_unstable,'k--','LineWidth',2)
plot(Om,sqrt(zna_max),'k--')
plot(Om,sqrt(zna_min),'k--')
 ax = gca;
ax.FontSize = 15; 
xlabel('$\Omega$','FontSize',16.5,'interpreter','latex') 
ylabel('$b$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis

