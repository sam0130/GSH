clc; clear; close all;
A = [0 0 1;
    0 0 0;
    1 0 0];

[vect, value] = eig(A);



B = [0 0 1;
    0 0 0;
    1 0 0];

R = vect;

EigenValue = R*B*R';

StrainRate = R'*EigenValue*R;
