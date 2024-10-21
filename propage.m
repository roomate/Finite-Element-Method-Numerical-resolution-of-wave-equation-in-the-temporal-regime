function [Us, Kinetic, Potential, Times] = propage(M, K, interpU0, interpU1, dt, niter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propage :
% Propage les conditions initiales � partir du propagateur discret saute-mouton.
%
% SYNOPSIS [Us, Kinetic, Potential, Times] = propage(M, K, interpU0, interpU1, dt, niter)
%
% INPUT  * M, K : matrice de masse et de rigidit�
%        * interpU0, interpU1 : interpol�e en espace des conditions initiales
%        * dt : pas de temps du sch�ma.
%        * niter : nombre d'it�rations
%
% OUTPUT * Us : saolution discr�te � tous les pas de temps
%               avec les conditions initiales (matrice Nbpt x niter + 2)
%        * Kinetic, Potential : �nergies cin�tique et potentielle (vecteur niter)
%        * Times : vecteur temps (vecteur niter)
%
% NOTE (1) la matrice M^{-1} K est calcul�e et stock�e avant la boucle en
% temps (peu efficace).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calcul et stockage de la matrice M^{-1} K.
% A COMPLETER

% Allocation.
Nbpts = length(interpU0);
Us = zeros(Nbpts, niter + 2);
Kinetic = zeros(niter, 1);
Potential = zeros(niter, 1);
Times = zeros(niter, 1);

% Conditions initiales.
Us(:, 1) = interpU0;
Us(:, 2) = 

% Iteriations.
for i = 1:niter
    
    % Calcul des �nergies.
    % A COMPLETER
    
    % Calcul de la solution � l'it�ration suivante.
    % A COMPLETER
    
    Times(i) = i * dt;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021