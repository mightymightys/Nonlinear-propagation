% ftAxis
% Création de deux tableaux (fréquence et temps) calibrés l'un par rapport à l'autre pour
% une FFT. L'espacement entre deux points de fréquence est l'inverse de la largeur totale
% dans l'espace des temps. Inversement, l'espacement entre deux points de temps est l'inverse
% de la largeur totale en fréquence.
% Les axes de fréquence et de temps sont supposés centré sur zéro
%
% nPoints : Nombre de points (de préférence une puissance de 2)
% nuMax   : Valeur maximale de la fréquence (ie moitié de la largeur totale)
% nu : Axe des fréquences. Le point d'indice nPoints/2 vaut toujours zéro
% t  : Axe des temps. Le point d'indice nPoints/2 vaut toujours zéro
function [nu, t] = ftAxis(nPoints, nuMax)
  deltaNu = 2*nuMax/nPoints;
  deltaT = 1/(2*nuMax);
  nu = -nuMax:2*nuMax/nPoints:nuMax-(2*nuMax/nPoints);
  t = -nPoints/2*deltaT:deltaT:(nPoints/2-1)*deltaT;
