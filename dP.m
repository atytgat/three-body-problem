function [out] = dP( a,x1i,x2i,x3i,r12,r13,r23 )
m1 = 1;                                  % masse du soleil [masse solaire]
m2 = 9.542986425 * 10^-4;                     % masse de jupiter [masse solaire]
m3 = 2.857315234 * 10^-4;                  % masse de Saturne [masse solaire] 
G =  0.01720209895^2;                    % constante grav de Gauss [AU^(-3/2) * jour solaire moyen * (masse solaire)^(1/2)]

if a == 1
    
    out = -G*m1*m2*(x1i - x2i) / r12^3 - G*m1*m3*(x1i - x3i) / r13^3; % Dérivée de l'impulsion du Soleil
    
elseif a == 2
    
    out = G*m1*m2*(x1i - x2i) / r12^3 - G*m2*m3*(x2i - x3i) / r23^3; % Dérivée de l'impulsion de Jupiter
    
else 
    
    out = G*m1*m3*(x1i - x3i) / r13^3 + G*m2*m3*(x2i - x3i) / r23^3; % Dérivée de l'impulsion de Saturne
    
end
end


