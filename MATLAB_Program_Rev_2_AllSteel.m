%% Introduction
%{
Created by Michael Stewart
Mar 27 2021
Student Number: 214693626

Revision 1

Notes about Running Program:



SELF CONTAINED VERSION
This version of the Program is self-contained. Meaning that it does not use
any external functions as requested per lab Manual. This will be changed in
future iterations. 

Reusibility: The calcualtions section of this code is made to be completely
reusable. The Intialization section (Setting the parameters of the
Elements)is not Reusable as the possibility in element lengths,
angles, interconnectedness, and Young's modulus is infinte. Also note that
the Visualization section of the program is only used to help the user to
visualize the final design chosen. 

The Approach taken to create this program was to hardcode Arrays of Element
and Node Parameters and create functions and other calculations that can
input these arrays to give proper answers

Assumptions:
 
Counter Clockwise Rotation is Positive
Clockwise Rotation is Negative
Units are S.I. (m, Pa, N, etc), and Angles are in Radians


%}
clear all

%% Visualization Section
%Visual Representaiton of Chosen Final Design (Non Reusable)

%This Section helps to visualize the created structure.


%Note that this part of the code is not Reusable and Cannot be used to
%Solve different problems

% X-Y Plane of the Sructure
x=linspace(0,3);y=-x+3; 
x2=linspace(0,2); y2=ones(1,length(x2));
x3=linspace(0,1);y4=-x3+1; y5=-x2+2;
figure(1);

plot(x,y,'r',x2,y2,'b',x3,2*y2,'b',2*y2,x3,'b',y2,x2,'b',x3,y4,'r',x2,y5,'r',zeros(1,length(x)),x,'b',x,zeros(1,length(x)),'r')
title('Visual Representation of Structure');
legend('Steel Beams','Proposed Aluminum Beams');
xlabel('X Axis (m)'); ylabel('Y Axis (m)');

% Cross Sectional Area
x=linspace(0,0.004); x1=linspace(0,0.008);
y=ones(1,length(x)); 
figure(2);
plot(x,0.008.*y,x,0.*y,0.*y,x1,0.004.*y,x1);
xlim([-0.01,0.01]); ylim([-0.01,0.01]);
title('Cross Sectional Area of the Beam (in meters)');
xlabel('Width(m)'); ylabel('Height (m)');
%% Intialization Section:
%Setting up Elemental Values (Non Reusable)

A=0.008*0.004; %m 



Alum_E=69e09; %Pascals
Steel_E=200e+09; %Pascals

% Element 1 (between Node 1 and 2
Length_1=1; %m
Angle_1= 0; %Degrees and Radians
E_1=Steel_E;
Node_connect_1=[1,2];

% Element 2(between Node 2 and 3)
Length_2=1; %m
Angle_2=0;
E_2=Steel_E;
Node_connect_2=[2,3];

%Element 3(between Node 3 and 4)
Length_3=1;
Angle_3=0;
E_3=Steel_E;
Node_connect_3=[3,4];


%Element 4 (between Node 1 and 5)
Length_4=1;
Angle_4=deg2rad(90); %Covnerts -180 degrees to radians
E_4=Steel_E;
Node_connect_4=[1,5];

%Element 5(between Node 2 and 5)
Length_5=sqrt(2); %comes from sqrt(1^2+1^2)
Angle_5=deg2rad(135); % Converts degs to radians
E_5=Steel_E;
Node_connect_5=[2,5];

%Element 6 (between Node 2 and 6)
Length_6=1;
Angle_6=deg2rad(90);
E_6=Steel_E;
Node_connect_6=[2,6];

%Element 7 (between 3 and 6)
Length_7=sqrt(2);
Angle_7=deg2rad(135);
E_7=Steel_E;
Node_connect_7=[3,6];


%Element 8 (between Node 2 and 7)
Length_8=1;
Angle_8=deg2rad(90);
E_8=Steel_E;
Node_connect_8=[3,7];

%Element 9(between 4 and 7)
Length_9=sqrt(2);
Angle_9=deg2rad(135);
E_9=Steel_E;
Node_connect_9=[4,7];

%Element 10(between 5 and 6)
Length_10=1;
Angle_10=0;
E_10=Steel_E;
Node_connect_10=[5,6];

%Element 11(between 6 and 7)
Length_11=1;
Angle_11=0;
E_11=Steel_E;
Node_connect_11=[6,7];

%Element 12(between 5 and 8)
Length_12=1;
Angle_12=deg2rad(90);
E_12=Steel_E;
Node_connect_12=[5,8];

%Element 13(between 6 and 8)
Length_13=sqrt(2);
Angle_13=deg2rad(135);
E_13=Steel_E;
Node_connect_13=[6,8];

%Element 14(between 6 and 9)
Length_14=1;
Angle_14=deg2rad(90);
E_14=Steel_E;
Node_connect_14=[6,9];

%Element 15(between 7 and 9)
Length_15=sqrt(2);
Angle_15=deg2rad(135);
E_15=Steel_E;
Node_connect_15=[7,9];


%Element 16(between 8 and 9)
Length_16=1;
Angle_16=deg2rad(0);
E_16=Steel_E;
Node_connect_16=[8,9];

%Element 17 (between 8 and 10)
Length_17=1;
Angle_17=deg2rad(90);
E_17=Steel_E;
Node_connect_17=[8,10];

%Element 18(Node betwen 9 to 10)
Length_18=sqrt(2);
Angle_18=deg2rad(135);
E_18=Steel_E;
Node_connect_18=[9,10];


%Hard Coded Arrays. To be used as an Example

Length_Array=[Length_1,Length_2,Length_3,Length_4,...
    Length_5,Length_6,Length_7,Length_8,Length_9,Length_10,...
    Length_11,Length_12,Length_13,Length_14,Length_15,Length_16,Length_17,Length_18];

E_Array=[E_1,E_2,E_3,E_4,E_5,E_6,E_7,E_8,E_9,E_10,E_11,E_12,E_13,E_14,E_15,E_16,E_17,E_18];

Angle_Array=[Angle_1,Angle_2,Angle_3,Angle_4,Angle_5,Angle_6,...
    Angle_7,Angle_8,Angle_9,Angle_10,Angle_11,Angle_12,Angle_13,Angle_14...
    ,Angle_15,Angle_16,Angle_17,Angle_18];

Num_of_elements=length(Length_Array);
A_Array=A.*ones(1,Num_of_elements); %Since all the Cross Sectional Areas are the same

Interconnect_mat=[Node_connect_1;Node_connect_2;Node_connect_3;...
    Node_connect_4;Node_connect_5;Node_connect_6;Node_connect_7;...
    Node_connect_8;Node_connect_9;Node_connect_10;Node_connect_11;...
    Node_connect_12;Node_connect_13;Node_connect_14;Node_connect_15;...
    Node_connect_16;Node_connect_17;Node_connect_18];





Node_mat=[0,0;1,0;2,0;3,0;0,1;1,1;2,1;0,2;1,2;0,3]; %Creates the Node Matrix

BC_Vector=[1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,1]; %Boundary Conditon Boolean Vector

GlobalForceVector=[0,0,0,-3000,0,-3000,0,-3000,0,0,0,0,0,0,0,0,0,0,0,0]; %Global Force Vector 

%% Generation of Global Stiffness Matrix

%This Section of the Program created the Global Stiffness Matrix

Zeros_mat=zeros((max(max(Interconnect_mat)))*2,(max(max(Interconnect_mat)))*2,length(Angle_Array));
%Takes the Maximum Node Value in the Interconnection Matrix and Multiplies
%it by 2

% In our Example the Zero Matrix is a 20 by 20 by 18
%This is because we need to represent X values, Y values in 10 Nodes with
%18 Elements


Global_Stiffness_2=Zeros_mat;

for i =1:length(Angle_Array) %iterations through the Elements 
    
    Global_Stiffness_1(:,:,i)=ElementStiffnessMatrixGlobal(E_Array(i),A_Array(i),Length_Array(i),Angle_Array(i));
    % Each 'Slice' is the Elements Stiffness Matrix
    
    %The Next Part is to Move the Values found in Global Stiffness_1 into a
    %big 2D array, one that could easily be added together to get the
    %global matrix
    
    
    %We must define our indexes of the X(i,1)-X(i,1) Values:
    Index_1=Interconnect_mat(i,1)*2-1; Index_2=Interconnect_mat(i,2)*2-1;
    %These index values use the attached nodel numbers.
    %in our example our Zeroes Matrix is 20 by 20 by 18 since we have 10
    %Nodes with X and Y values for the Nodes and 18 Elements. The Index
    %values use the Node's number multiplies it by 2 and subtracts 1 to
    %represent the Global Stiffness Matrix's Index Value in 2D
    
        
    %We now assign the Global Stiffness Matrix values to the proper indices
    
    %Qudrant 2 (Top Left):
    
    Global_Stiffness_2(Index_1,Index_1,i)=Global_Stiffness_1(1,1,i);
    %X1-X1 Component
    Global_Stiffness_2(Index_1+1,Index_1,i)=Global_Stiffness_1(2,1,i);
    %X1-Y1 Component
    Global_Stiffness_2(Index_1,Index_1+1,i)=Global_Stiffness_1(1,2,i);
    %Y1-X1 Component
    Global_Stiffness_2(Index_1+1,Index_1+1,i)=Global_Stiffness_1(2,2,i);
    %Y1-Y1 Component
    
    %Quadrant 4 (Bottom Right)
    Global_Stiffness_2(Index_2,Index_2,i)=Global_Stiffness_1(3,3,i);
    %X2-X2 Component
    Global_Stiffness_2(Index_2+1,Index_2,i)=Global_Stiffness_1(4,3,i);
    %X2-Y2 Component
    Global_Stiffness_2(Index_2,Index_2+1,i)=Global_Stiffness_1(3,4,i);
    %Y2-X2 Component
    Global_Stiffness_2(Index_2+1,Index_2+1,i)=Global_Stiffness_1(4,4,i);
    %Y2-Y2 Component
    
    %Quadrant 3 (Bottom Left)
    Global_Stiffness_2(Index_1,Index_2,i)=Global_Stiffness_1(1,3,i);
    %X1-X2 Component
    Global_Stiffness_2(Index_1+1,Index_2,i)=Global_Stiffness_1(2,3,i);
    %X1-Y2 Component
    Global_Stiffness_2(Index_1,Index_2+1,i)=Global_Stiffness_1(1,4,i);
    %Y1-X2 Component
    Global_Stiffness_2(Index_1+1,Index_2+1,i)=Global_Stiffness_1(2,4,i);
    %Y1-Y2 Component
    
    %Quadrant 1 (Top Right)
    Global_Stiffness_2(Index_2,Index_1,i)=Global_Stiffness_1(3,1,i);
    %X2-X1 Component
    Global_Stiffness_2(Index_2,Index_1+1,i)=Global_Stiffness_1(3,2,i);
    %X2-Y1 Component
    Global_Stiffness_2(Index_2+1,Index_1,i)=Global_Stiffness_1(4,1,i);
    %Y2-X1 Component
    Global_Stiffness_2(Index_2+1,Index_1+1,i)=Global_Stiffness_1(4,2,i);
    %Y2-Y1 Component   
end
%Now that we have Propagted our Global Stiffness Matrices for each Element
%, we can add up all the 'Slices' to get a 2D Sum-- This is our Global
%Stiffness Matrix!
GlobalStiffnessMatrix_Example=sum(Global_Stiffness_2,3);




%% Applying Boundary Conditions

%Applying BC's

%Applying BC to Global Stiffness Matrix
List_Non_zeros=find(BC_Vector); % Finds the indexes of the Boundary Conditions
%finds the indexs of Boundary Conditons

BC_GSM=GlobalStiffnessMatrix_Example;

BC_GSM([List_Non_zeros],:)=[];
BC_GSM(:,[List_Non_zeros])=[];
%deletes the rows and columns that have Boundary Conditions

%Applying BC to Forces Vector
BC_Forces=GlobalForceVector;
BC_Forces([List_Non_zeros])=[];
%Deletes Froces Values which have bounary conditions

%% Solve for Displacements

Displacements=inv(BC_GSM)*(transpose(BC_Forces));
%This gives the Displacements in meters

%% Solve for Reaction Forces
List_zeros=find(~BC_Vector); %Finds index of Zeros in BC Unpacked Array


RF_GSM=GlobalStiffnessMatrix_Example;
RF_GSM(:,[List_Non_zeros])=[]; %Deletes Columns with ones (that have boundary conditons)
RF_GSM([List_zeros],:)=[]; %Deletes Rows with Zeros (that are not boundary conditioned)

Reaction_force=(RF_GSM)*(Displacements);


%% Create Full Arrays of the Reaction Forces and Displacements for each Node

%This Section populates the Reaction and Displacement array so that they
%have a value of 0 when they are not considered in the above calculations
index=1;
index2=1;

for i=1:length(BC_Vector)
    if (BC_Vector(i)==1)
        Total_displacements(i)=0;
        ReactionForce_total(i)=Reaction_force(index2);
        index2=index2+1;
    else
        Total_displacements(i)=Displacements(index);
        index=index+1;
    end
end


%% Calculate Stress of Each Element
for i=1:length(Length_Array)
    
    StressAngleArray=[-cos(Angle_Array(i)),-sin(Angle_Array(i)),...
        cos(Angle_Array(i)),sin(Angle_Array(i))];
    X_1=Interconnect_mat(i,1)*2-1; X_2=Interconnect_mat(i,2)*2-1; %indexs
    Y_1=Interconnect_mat(i,1)*2; Y_2=Interconnect_mat(i,2)*2; 
    %Using the Node Numbers mentioned in the interconnectiveness matrix, we
    %come up with indexes to be accessed for each element
    
    Element_displacement=[Total_displacements(X_1);Total_displacements(Y_1);...
        Total_displacements(X_2);Total_displacements(Y_2)];
    Stress(i)=(StressAngleArray*Element_displacement)*E_Array(i)/Length_Array(i); 
    %Produces Stresses of the Elements
    %Using the Angle Array we can multiply it by the displacements
end

%% Output/Display Section
%Plot Deformed and Undeformed Nodes


%Repack Displacements in x and y coordinates
Displacements_Packed(:,1)=Total_displacements(1:2:end);
Displacements_Packed(:,2)=Total_displacements(2:2:end);
%this is done because it's easier to add to the exisitng Node matrix

Scale=20;

Displaced_Nodes=(Node_mat+Scale.*(Displacements_Packed));
%Adds the existing x,y corrdinates of the element to the node displacements

figure(3)
plot(1:2)=scatter(Node_mat(:,1),Node_mat(:,2),'filled','b');hold on
plot(3:4)=scatter(Displaced_Nodes(:,1),Displaced_Nodes(:,2),'r','filled');
hold off;
legend(plot([2,4]) ,{'Orginal Nodes','Displaced Nodes (Scale 20)'})
xlabel('X Axis'); ylabel('Y Axis')

% Print Deformation and Stress
disp('Reaction Forces(N):')
disp (transpose(ReactionForce_total))
disp('The Element Stress is (Pa):')
disp(transpose(Stress))
disp('The Nodal Displacements are (m)[x,y]:')
disp(Displacements_Packed)

%% Local Function
function [output]= ElementStiffnessMatrixGlobal(E,A,L,Angle)

%First Step is to Calculate the Elemental Stiffness Matrix

%Matrix for Elemental Stiffness Matrix
mat=[1,0,-1,0;0,0,0,0;-1,0,1,0;0,0,0,0];
E_Stiffness_mat=mat*E*A/L;
%Second Step is to Calculate the Transformation Matrix of Each Elemnet
Tranform_Matrix=[cos(Angle),sin(Angle),0,0;sin(Angle),cos(Angle),0,0;...
     0,0,cos(Angle),sin(Angle);0,0,sin(Angle),cos(Angle)];
    %Each 'Slice' Represents the Element's Transformation Matrix
    
    %This Global Stiffness Matrix represents the Values of Each Element
    %which have not been added together.
output=transpose(Tranform_Matrix)*...
        E_Stiffness_mat*Tranform_Matrix;


end


%%% End of Program %%%  