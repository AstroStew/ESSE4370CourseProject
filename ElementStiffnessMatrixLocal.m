function [output]= ElementStiffnessMatrixLocal(E,A,L)
%Matrix for Elemental Stiffness Matrix
mat=[1,0,-1,0;0,0,0,0;-1,0,1,0;0,0,0,0];
output=mat*E*A/L;
end
