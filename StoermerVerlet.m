A = 365*500;            % interval de temps
k = 10;                  % pas de temps/précision
N = 1 + round(A / k);    % nbr de points de discrétisation

m1 = 1;                                  % masse du soleil [masse solaire]
m2 = 9.542986425 * 10^-4;                     % masse de jupiter [masse solaire]
m3 = 2.857315234 * 10^-4;                  % masse de Saturne [masse solaire] 
G =  0.01720209895^2;                    % constante grav de Gauss [AU^(-3/2) * jour solaire moyen * (masse solaire)^(1/2)]

% Le 2017-12-01T20:55:29.00
% Les positions sont exprimées en UA et les impulsions en UA * masse solaire /jour

% CI Soleil
x10 = 0;                
y10 = 0;
z10 = 0;
px10 = -0.0043508739213 * m2 - 0.0052777229823 * m3*0;          
py10 = 0.0052486922392 * m2 - 0.0000052313702 * m3*0;  
pz10 = 0.0023556626834 * m2 + 0.0002250824195 * m3*0;

% CI Jupiter
x20 = -4.3960041666534;     
y20 = -2.9773474501052;
z20 = -1.1691562160985;
px20 = 0.0043508739213 * m2;          
py20 = -0.0052486922392 * m2;  
pz20 = -0.0023556626834 * m2;

% CI Saturne
x30 = -0.1129204020027;     
y30 = -9.3029902925425;
z30 = -3.8378234683606;
px30 = 0.0052777229823 * m3;          
py30 = 0.0000052313702 * m3;  
pz30 = -0.0002250824195 * m3;

% Soleil
x1 = zeros(N,3);
x1(1,1) = x10;
x1(1,2) = y10;
x1(1,3) = z10;

p1 = zeros(N,3);
p1(1,1) = px10;
p1(1,2) = py10;
p1(1,3) = pz10;

% Jupiter
x2 = zeros(N,3);
x2(1,1) = x20;
x2(1,2) = y20;
x2(1,3) = z20;

p2 = zeros(N,3);
p2(1,1) = px20;
p2(1,2) = py20;
p2(1,3) = pz20;

% Saturne 
x3 = zeros(N,3);
x3(1,1) = x30;
x3(1,2) = y30;
x3(1,3) = z30;

p3 = zeros(N,3);
p3(1,1) = px30;
p3(1,2) = py30;
p3(1,3) = pz30;

% Prédictions
p1pred = zeros(1,3);
p2pred = zeros(1,3);
p3pred = zeros(1,3);

% Grandeurs caractéristiques du système
Energie = zeros(N,1);
Impulsion = zeros(N,1);
MomentAngulaire = zeros(N,1);

for s = 1:N
    
    % Distances entre les corps
    r12 = norm(x1(s,:)-x2(s,:)); % Soleil et Jupiter
    r13 = norm(x1(s,:)-x3(s,:)); % Soleil et Saturne
    r23 = norm(x2(s,:)-x3(s,:)); % Jupiter et Saturne
    
    for i=1:3
        
        % Prédictions des impulsions
        p1pred(1,i) = p1(s,i) + k/2*dP(1,x1(s,i),x2(s,i),x3(s,i),r12,r13,r23);
        p2pred(1,i) = p2(s,i) + k/2*dP(2,x1(s,i),x2(s,i),x3(s,i),r12,r13,r23);
        p3pred(1,i) = p3(s,i) + k/2*dP(3,x1(s,i),x2(s,i),x3(s,i),r12,r13,r23);
        
        % Nouvelles positions
        x1(s+1,i) = x1(s,i) + k*dQ(p1pred(1,i),m1);
        x2(s+1,i) = x2(s,i) + k*dQ(p2pred(1,i),m2);
        x3(s+1,i) = x3(s,i) + k*dQ(p3pred(1,i),m3);
        
    end
    
    % Nouvelles distances
    r12pred = norm(x1(s+1,:)-x2(s+1,:)); % Soleil et Jupiter
    r13pred = norm(x1(s+1,:)-x3(s+1,:)); % Soleil et Saturne
    r23pred = norm(x2(s+1,:)-x3(s+1,:)); % Jupiter et Saturne

    for i=1:3
        
        % Nouvelles impulsions
        p1(s+1,i) = p1pred(1,i) + k/2*dP(1,x1(s+1,i),x2(s+1,i),x3(s+1,i),r12pred,r13pred,r23pred);
        p2(s+1,i) = p2pred(1,i) + k/2*dP(2,x1(s+1,i),x2(s+1,i),x3(s+1,i),r12pred,r13pred,r23pred);
        p3(s+1,i) = p3pred(1,i) + k/2*dP(3,x1(s+1,i),x2(s+1,i),x3(s+1,i),r12pred,r13pred,r23pred);
        
    end



    Impulsion(s) = norm(p1(s,:)) + norm(p2(s,:)) + norm(p3(s,:));
    MomentAngulaire(s) = norm(x1(s,:)) * norm(p1(s,:)) + norm(x2(s,:)) * norm(p2(s,:)) + norm(x3(s,:)) * norm(p3(s,:));
    Energie(s) = norm(p1(s,:))^2 / (2*m1) + norm(p2(s,:))^2 / (2*m2) + norm(p3(s,:))^2 / (2*m3) - G*m1*m2/r12 - G*m1*m3/r13 - G*m2*m3/r23;
    
end


figure(1)
plot3(x1(:,1),x1(:,2),x1(:,3),'r',x2(:,1),x2(:,2),x2(:,3),'b',x3(:,1),x3(:,2),x3(:,3),'g')
title('Schéma de Stoermer-Verlet')
xlabel('x')
ylabel('y')
zlabel('z')
legend('Soleil', 'Jupiter', 'Saturne')

figure(2)
plot(1:N,Energie)
title('Energie totale')
xlabel('Temps')
ylabel('Energie')

figure(3)
plot(1:N,Impulsion,'r')
plot(1:N,Impulsion,'r')
title('Impulsion totale')
xlabel('Temps')
ylabel('Norme de l impulsion')

figure(4)
plot(1:N,MomentAngulaire)
title('Moment angulaire total')
xlabel('Temps')
ylabel('Norme du moment angulaire')
