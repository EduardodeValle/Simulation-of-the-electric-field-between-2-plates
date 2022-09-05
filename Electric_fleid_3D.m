clc
clear all
close all

%Valores de la placa
h1 = 1; %Altura de la placa 1
h2 = 1; %Altura de la placa 2
Qp = 1; %Carga de la placa 1 (positiva)
Qn = 1; %Carga de la placa 2 (negativa)
d = 1; %Distancia entre placas
p = 1; %Profundidad de la placa
w = 1; %Ancho de la placa
xCp = -d/2; %Centro de la placa positiva en x
yCp = 0; %Centro de la placa positiva en y
zCp = 0; %Centro de la placa positiva en z
xCn = d/2; %Centro de la placa negativa en x
yCn = 0; %Centro de la placa negativa en y
zCn = 0; %Centro de la placa negativa en z

%Matriz del campo electrico
N=11; %Valor de N
minX=xCp-w; maxX=xCn+w; %Limites en x
minY=yCp-p; maxY=yCn+p; %Limites en Y
minZ=zCn-2*h1; maxZ=zCn+2*h1; %Limites en Z

xG=linspace(minX,maxX,N); %Espaciamiento entre x
yG=linspace(minY,maxY,N); %Espaciamiento entre y
zG=linspace(minZ,maxZ,N); %Espaciamiento entre z
[X,Y,Z] = meshgrid(xG,yG,zG); %Grid de los vectores

%Obtencion del campo eléctrico
Ex = 0; %Campo inicial en x
Ey = 0; %Campo inicial en y
Ez = 0; %Campo inicial en z

[Ex,Ey,Ez] = campo(X,Y,Z,xCp,yCp,zCp,h1,Qp,Ex,Ey,Ez); %Calculo del campo de la placa 1
[Ex,Ey,Ez] = campo(X,Y,Z,xCn,yCn,zCn,h2,Qn,Ex,Ey,Ez); %Calculo del campo de la placa 2

quiver3(X,Y,Z,Ex,Ey,Ez,"r");

%Graficacion del campo vectorial
hold on
axis ([-1 1 -1 1 -1 1])
% Create UIAxes
title('Campo vectorial')
xlabel('X')
ylabel('Y')
zlabel('Z')

%Graficacion de las lineas con carga
for i = 0:.1:1
plot3(xCp,yCp,(zCp-h1*(i)),"o","Color","b"); %Graficacion de puntos de la placa
plot3(xCp,yCp,(zCp+h1*(i)),"o","Color","b");
plot3(xCn,yCn,(zCn-h2*(i)),"o","Color","b");
plot3(xCn,yCn,(zCn+h2*(i)),"o","Color","b");
end

posX = 0; %Obtener el valor de la coordenada en x
posY = 0; %Obtener el valor de la coordenada en y
posZ = 0.6; %Obtener el valor de la coordenada en z

Cx=0; Cy=0; Cz=0; %Campo inicial en x,y,z

[Cx,Cy,Cz] = campo_cords(posX,posY,posZ,xCp,yCp,zCp,h1,Qp,Cx,Cy,Cz); %Calculo del campo de la primera placa
[Cx,Cy,Cz] = campo_cords(posX,posY,posZ,xCn,yCn,zCn,h2,Qn,Cx,Cy,Cz); %Calculo del campo de la segunda placa

plot3(posX,posY,posZ,"o","Color","g");

function [Ex,Ey,Ez] = campo(X,Y,Z,xC,yC,zC,h,Q,Ex,Ey,Ez) %Calculo del campo vectorial
L=Q/h; %Lambda
eps=8.8541878176*10^(-12); %Constante de permitividad
kC = 1/(4*pi*eps); %Constante eléctrica
Rx = X - xC; %Desplazamiento en funcion del centro de la placa
Ry = Y - yC;
Rz = Z - zC;
sa2 = (h/2 - Rz) ./ (sqrt(Rx.^2 + Ry.^2 + (Rz - h/2).^2 )); %seno a2
sa1 = -(h/2 + Rz) ./ (sqrt(Rx.^2 + Ry.^2 + (Rz + h/2).^2 )); %seno a1
ca2 = (sqrt(Rx.^2 + Ry.^2)) ./ (sqrt(Rx.^2 + Ry.^2 + (Rz - h/2).^2 )); %coseno a2
ca1 = (sqrt(Rx.^2 + Ry.^2)) ./ (sqrt(Rx.^2 + Ry.^2 + (Rz + h/2).^2 )); %coseno a2
Ex = Ex +(kC .* ((L.*Rx) ./ (Rx.^2 + Ry.^2)) .* (sa2 - sa1)); %Campo en x
Ey = Ey + (kC .* ((L.*Ry) ./ (Rx.^2 + Ry.^2)) .* (sa2 - sa1)); %Campo en y
Ez = Ez + (kC .* ((L) ./ (sqrt(Rx.^2 + Ry.^2))) .* (ca2 - ca1)); %Campo en z
end

function [Ex,Ey,Ez] = campo_cords(X,Y,Z,xC,yC,zC,h,Q,Ex,Ey,Ez)
L=Q/h;
eps=8.8541878176*10^(-12);
kC = 1/(4*pi*eps);
Rx = X - xC;
Ry = Y - yC;
Rz = Z - zC;
sa2 = (h/2 - Rz) / (sqrt(Rx^2 + Ry^2 + (Rz - h/2)^2 ));
sa1 = -(h/2 + Rz) / (sqrt(Rx^2 + Ry^2 + (Rz + h/2)^2 ));
ca2 = (sqrt(Rx^2 + Ry^2)) / (sqrt(Rx^2 + Ry^2 + (Rz - h/2)^2 ));
ca1 = (sqrt(Rx^2 + Ry^2)) / (sqrt(Rx^2 + Ry^2 + (Rz + h/2)^2 ));
Ex = Ex +(kC * ((L*Rx) / (Rx^2 + Ry^2)) * (sa2 - sa1));
Ey = Ey + (kC * ((L*Ry) / (Rx^2 + Ry^2)) * (sa2 - sa1));
Ez = Ez + (kC * ((L) / (sqrt(Rx^2 + Ry^2))) * (ca2 - ca1));
end  