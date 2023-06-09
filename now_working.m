clear all
close all
clc

syms q8(t) q7(t) q6(t) q5(t) q4(t) q3(t) q2(t) q1(t) t a4 a2 a3 a1 qp1 qp2 qp3 qp4 qp5 qp6 qp7 qp8

RP=[0 0 0 0 0 0 0 0];

Q= [q1 q2 q3 q4 q5 q6 q7 q8];

Qp=[qp1 qp2 qp3 qp4 qp5 qp6 qp7 qp8];

GDL= size(RP,2);
GDL_str= num2str(GDL);

P(:,:,1)= [a1;0;0];
R(:,:,1)= [1 0 0;
           0 cos(q1) -sin(q1);
           0 sin(q1) cos(q1)]; 

P(:,:,2)= [0;0;0];
R(:,:,2)= [cos(q2) 0 sin(q2);
           0 1 0;
           -sin(q2) 0 cos(q2)];

P(:,:,3)= [0;0;0];
R(:,:,3)= [cos(q3) 0 sin(q3);
           0 1 0;
           -sin(q3) 0 cos(q3)];

P(:,:,4)= [0;0;0];
R(:,:,4)= [cos(q4) -sin(q4) 0;
          sin(q4) cos(q4) 0;
          0 0 1];

P(:,:,5)=[a2;0;0];
R(:,:,5)=[1 0 0;
          0 1 0;
          0 0 1];

P(:,:,6)=[0;0;0];
R(:,:,6)=[cos(q5) -sin(q5) 0;
          sin(q5) cos(q5) 0;
          0 0 1];

P(:,:,7)=[0;0;0];
R(:,:,7)=[1 0 0;
           0 cos(q6) -sin(q6);
           0 sin(q6) cos(q6)];

P(:,:,8)=[0;0;0];
R(:,:,8)=[cos(q7) 0 sin(q7);
            0 1 0;
          -sin(q7) 0 cos(q7)];

P(:,:,9)=[0;0;a4];
R(:,:,9)=[1 0 0;
          0 1 0;
          0 0 1];

Vector_Zeros= zeros(1, 3);

%Inicializamos las matrices de transformación Homogénea locales
A(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
%Inicializamos las matrices de transformación Homogénea globales
T(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
%Inicializamos las posiciones vistas desde el marco de referencia inercial
PO(:,:,GDL)= P(:,:,GDL); 
%Inicializamos las matrices de rotación vistas desde el marco de referencia inercial
RO(:,:,GDL)= R(:,:,GDL); 


for i = 1:GDL
    i_str= num2str(i);
   %disp(strcat('Matriz de Transformación local A', i_str));
    A(:,:,i)=simplify([R(:,:,i) P(:,:,i); Vector_Zeros 1]);
   %pretty (A(:,:,i));

   %Globales
    try
       T(:,:,i)= T(:,:,i-1)*A(:,:,i);
    catch
       T(:,:,i)= A(:,:,i);
    end
    disp(strcat('Matriz de Transformación global T', i_str));
    T(:,:,i)= simplify(T(:,:,i));
    pretty(T(:,:,i))

    RO(:,:,i)= T(1:3,1:3,i);
    PO(:,:,i)= T(1:3,4,i);
    %pretty(RO(:,:,i));
    %pretty(PO(:,:,i));
end

%Calculamos el jacobiano lineal de forma analítica
Jv_a(:,GDL)=PO(:,:,GDL);
Jw_a(:,GDL)=PO(:,:,GDL);

for k= 1:GDL
    if RP(k)==0 
       %Para las juntas de revolución
        try
            Jv_a(:,k)= cross(RO(:,3,k-1), PO(:,:,GDL)-PO(:,:,k-1));
            Jw_a(:,k)= RO(:,3,k-1);
        catch
            Jv_a(:,k)= cross([0,0,1], PO(:,:,GDL));%Matriz de rotación de 0 con respecto a 0 es la Matriz Identidad, la posición previa tambien será 0
            Jw_a(:,k)=[0,0,1];%Si no hay matriz de rotación previa se obtiene la Matriz identidad
         end
     else
%         %Para las juntas prismáticas
        try
            Jv_a(:,k)= RO(:,3,k-1);
        catch
            Jv_a(:,k)=[0,0,1];
        end
            Jw_a(:,k)=[0,0,0];
     end
 end    

Jv_a= simplify (Jv_a);
Jw_a= simplify (Jw_a);

disp('Velocidad lineal obtenida mediante el Jacobiano lineal');
V=simplify (Jv_a*Qp');
pretty(V);
disp('Velocidad angular obtenida mediante el Jacobiano angular');
W=simplify (Jw_a*Qp');
    pretty(W);