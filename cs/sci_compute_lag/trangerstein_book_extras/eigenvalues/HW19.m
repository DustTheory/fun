%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 29 October 2013
% Lauren Lowman
%
% HW #19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

tic

%set matrix A
D = importdata('psmigr_1.mtx');
m = D(1,1);
n = D(1,2);
nonzero = D(1,3);

A = zeros(m,n);

for k = 2:nonzero+1
    i = D(k,1);
    j = D(k,2);
    A(i,j) = D(k,3);
end

toc

clear D;

% A = [2,2;3,5];
% [m,n] = size(A);

%compute SVD of A
%[U,S,V] = svd(A,'econ');
tic
[U,S,V] = svd(A);
toc

%evaluate norms
tic
P1 = U'*U - eye(m);
FN1 = norm(P1,'fro');
disp('FNORM1 = ')
disp(FN1)

%evaluate norms
P2 = V'*V - eye(m);
FN2 = norm(P2,'fro');
disp('FNORM2 = ')
disp(FN2)

%evaluate norms
P3 = A - U*S*V';
FN3 = norm(P3,'fro');
disp('FNORM3 = ')
disp(FN3)

toc
