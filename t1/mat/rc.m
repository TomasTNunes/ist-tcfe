close all
clear all

%%%%%%%%%%%%%%%%%%%EXAMPLE SYMBOLIC COMPUTATIONS%%%%%%%%%%%%%%%%%%%

%pkg load symbolic

R1 = 1.02943118797 
R2 = 2.01395929215 
R3 = 3.04865352258 
R4 = 4.00813564187 
R5 = 3.14877356154 
R6 = 2.0324721756 
R7 = 1.02032682556 
Va = 5.05481864136 
Id = 1.02872173547 
Kb = 7.12052712169 
Kc = 8.12923323408 

G1 = 1/R1
G2 = 1/R2
G3 = 1/R3
G4 = 1/R4
G5 = 1/R5
G6 = 1/R6
G7 = 1/R7


%%%%%Mesh method%%%%%

A = [R1+R3+R4,R3,R4,0;R4,0,R4+R6+R7-Kc,0;Kb*R3,Kb*R3-1,0,0;0,0,0,1]
B = [Va,0,0,Id]'
C = A\B

Ise = C(1)
Isd = C(2)
Iie = C(3)
Iid = C(4)

V1m = Va
V2m = V1m-R1*Ise
V3m = V2m+R2*Isd
V4m = R4*(Ise+Iie)
V5m = V4m-R5*(Isd-Iid)
V6m = -R6*Iie
V7m = V6m-R7*Iie
V8m = V6m


%%%%%Nodal method%%%%%

D = [-G1,G1+G2+G3,-G2,-G3,0,0,0;G1,-G1,0,-G4,0,-G6,0;0,Kb,0,-Kb-G5,G5,0,0;G1,-G1,0,-G4,0,G7,-G7;1,0,0,0,0,0,0;0,0,0,1,0,Kc*G6,-1;0,Kb+G2,-G2,-Kb,0,0,0]
E = [0,0,Id,0,Va,0,0]'
F = D\E

V1n = F(1)
V2n = F(2)
V3n = F(3)
V4n = F(4)
V5n = F(5)
V6n = F(6)
V7n = F(7)
V8n = V6n



