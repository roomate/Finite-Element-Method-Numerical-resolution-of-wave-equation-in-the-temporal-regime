function [Us, Kinetic, Potential, Times] = propage_cond(M, K, interpU0, interpU1, dt, niter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propage_cond :
% Propage les conditions initiales � partir du propagateur discret saute-mouton.
%
% SYNOPSIS [Us, Kinetic, Potential, Times] = propage_cond(M, K, interpU0, interpU1, dt, niter)
%
% INPUT  * M, K : matrice de masse condens�e et de rigidit�
%        * interpU0, interpU1 : interpol�e en espace des conditions initiales
%        * dt : pas de temps du sch�ma.
%        * niter : nombre d'it�rations
%
% OUTPUT * Us : solution discr�te � tous les pas de temps
%               avec les conditions initiales (matrice Nbpt x (niter + 2))
%        * Kinetic, Potential : �nergies cin�tique et potentielle (vecteur niter)
%        * Times : vecteur temps (vecteur niter)
%
% NOTE (1) On suppose que la matrice de masse est diagonale de sorte que
% l'op�ration M \ b est efficace.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocation.
Nbpt = length(interpU0);
Us = zeros(Nbpt, niter + 2);
Kinetic = zeros(niter, 1);
Potential = zeros(niter, 1);
Times = zeros(niter, 1);
Z = M\K;
% Conditions initiales.
Us(:, 1) = interpU0;
Us(:, 2) = (eye(Nbpt) - (dt^2)/2*Z)*interpU0 + dt*interpU1;


% Interiations.
for i = 1:niter
    
    % Calcul des �nergies.
    deltaV = Us(:,i+1) - Us(:,i);
    Kinetic(i) = 1/2*(deltaV)'/dt*(M - K*(dt^2 / 4))*(deltaV)/dt;
    Potential(i) = 1/8*(Us(:,i+1) + Us(:,i))'*K*(Us(:,i+1) + Us(:,i));
    % Calcul de la solution par r�solution directe du syst�me lin�aire
    % (suppose que l'op�ration M \ b soit efficace, i.e. M diagonale).
    Us(:,i+2) = (2*eye(Nbpt,'like',K) - (dt^2)*Z)*Us(:,i+1) - Us(:,i) ;
    
    Times(i) = i * dt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021