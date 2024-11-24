tic
clear all;
close all;
clc

% AJUSTES RAPIDOS
% Los valores mas pequeños tienen mayor efecto

% VELOCIDAD DE METEORITOS
vel = 2;
% INTENSIDAD DE METEORITOS
intens = 5;

%Se carga laberinto de la imagen Creación del ambiente
lab = imread('mapa.png');
labB = 1 > lab;
labB = labB(:,:,1)
labBi = zeros(100,100)
labBi(labB)=255
%imtool(labBi);

%se instancian las estaciones
labBi(20,44)=1;
labBi(41,48)=1;
labBi(63,48)=1;
labBi(76,44)=1;

D = labBi;

initialPos = [[20,44];[41,48];[63,48];[76,44]];

NA=4

M=zeros(1,NA);
N=zeros(1,NA);

A.pos=zeros(size(labBi));
A.status=zeros(size(labBi)); %%Atacando, regresando y esperando
A.meteoritos=zeros(size(labBi));

%inicializacion de los agentes
for k=1:NA
    
    x=initialPos(k,2);
    y=initialPos(k,1);
    
    A.pos(y,x)=1;
    A.status(y,x)=0;
    
    M(k)=y; 
    N(k)=x;
    
end

%Numero de iteraciones 
it=2500;

%3. Se inician las iteraciones y se escriben las reglas
%A.status 0 indica espera
%A.status 1 indica ataque
%A.status 2 indica regresando

for itera=1:it
    if(mod(itera,intens) == 0)
        row = randi([21,84]);
        col = 99;
        A.meteoritos(row,col) = 1000;
    end

    % Movimiento de los meteoritos
    if (mod(itera,vel)==0)
        meteoritos = find(A.meteoritos == 1000);  % Encuentra las posiciones de los meteoritos
        if ~isempty(meteoritos)
            for i = 1:length(meteoritos)
                [meteorito_y, meteorito_x] = ind2sub(size(D), meteoritos(i));

                % Mover meteorito hacia la izquierda
                if meteorito_x > 2
                    A.meteoritos(meteorito_y, meteorito_x) = 0;  % Limpiar la posición actual
                    if(D(meteorito_y, meteorito_x - 1) == 255)
                        D(meteorito_y-1:meteorito_y+1,meteorito_x-1:meteorito_x+1) = 0;
                        labBi(meteorito_y-1:meteorito_y+1,meteorito_x-1:meteorito_x+1) = 0;
                    else
                        A.meteoritos(meteorito_y, meteorito_x - 1) = 1000;  % Mover meteorito una columna a la izquierda
                    end
                else
                    A.meteoritos(meteorito_y, meteorito_x) = 0;  % Eliminar meteorito al llegar al borde izquierdo
                end
            end
        end
    end

    for h=1:length(M)
        %almacenamiento de las coordenadas originales del agente al inicio de cada iteración.   
        posm=M(h);
        posn=N(h);

        Block=A.meteoritos(M(h)-10:M(h)+10,N(h)-10:N(h)+10);
        Lo=ismember(Block,1000);

        % Valida que el agente no tiene un meteorito en su vecindario, pero
        % solo se activa si el agente esta en espera
        if (sum(Lo(:))>0 && A.status(M(h),N(h))==0)
            % si hay un meteorito, la nave pasa a atacar
            A.status(M(h),N(h))=1;
        end

        % si se esta regresando
        if (A.status(posm,posn)==2)
            %Se deja el espacio anterior vacio
            A.pos(posm,posn)=0;
            A.status(posm,posn)=0;

            if (initialPos(h,2) > posn)
                posn = posn + 1;
            elseif (initialPos(h,2) < posn)
                posn = posn - 1;
            end

            if (initialPos(h,1) > posm)
                posm = posm + 1;
            elseif (initialPos(h,1) < posm)
                posm = posm - 1;
            end

            if(D(posm,posn) == 1)
                A.status(posm,posn) = 0;
            else

                A.status(posm,posn)=2;
            end

            A.pos(posm,posn)=1;

        end
                
        % Si esta en ataque
        if (A.status(posm,posn)==1)
            [meteorito_pos_y,meteorito_pos_x] = find(Block == 1000);
            if([meteorito_pos_y,meteorito_pos_x])
                meteorito_pos_x = meteorito_pos_x(1);
                meteorito_pos_y = meteorito_pos_y(1);

                A.pos(posm,posn)=0;
                A.status(posm,posn)=0;

                if (meteorito_pos_x > 11)
                    posn = posn + 1;
                elseif (meteorito_pos_x < 11)
                    posn = posn - 1;
                end

                if (meteorito_pos_y > 11)
                    posm = posm + 1;
                elseif (meteorito_pos_y < 11)
                    posm = posm - 1;
                end

                if(A.meteoritos(posm,posn) == 1000)
                    A.meteoritos(posm,posn) = 0;
                    A.status(posm,posn)=2;
                else
                    A.status(posm,posn)=1;
                end

                A.pos(posm,posn)=1;
            else
                A.status(posm,posn)=2;
            end
        end

        M(h) = posm;
        N(h) = posn;
    end

end

toc