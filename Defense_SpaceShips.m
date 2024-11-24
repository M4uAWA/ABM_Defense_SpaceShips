tic
clear all;
close all;
clc

%Se carga laberinto de la imagen Creación del ambiente
lab = imread('mapa.png');
labB = 1 > lab;
labB = labB(:,:,1)
labBi = zeros(100,100)
labBi(labB)=255
%imtool(labBi);

%se inicializa el mapeo del area
Ib=zeros(100);

%se instancian las naves
Ib(20,44)=1;
Ib(41,48)=1;
Ib(63,48)=1;
Ib(76,44)=1;

initialPos = [[20,44];[41,48];[63,48];[76,44]];

Chamfer; %mapa de distancias, recibe la imagen y regresa distancias euclidianas o con el algoritmo de manhattan
%Se añade laberinto a chamfer

D(labBi==255)=255;

[m, n]=size(Ib)

NA=4

M=zeros(1,NA);
N=zeros(1,NA);

A.pos=zeros(size(Ib));
A.status=zeros(size(Ib)); %%Atacando, regresando y esperando
A.meteoritos=zeros(size(Ib));

%meteorito fijo de prueba

A.meteoritos(27,50) = 1000;
A.meteoritos(48,50) = 1000;
A.meteoritos(28,49) = 1000;
A.meteoritos(68,47) = 1000;

for k=1:NA
  
    x=initialPos(k,2);
    y=initialPos(k,1);
    
    A.pos(y,x)=1;
    A.status(y,x)=0;

    M(k)=y; 
    N(k)=x;
    
end

%Heatmap
num_rows = 100;
num_cols = 100;

% Crear el colormap de difuminado de colores: rojo, naranja, amarillo, verde, blanco, negro
cmap = [
    0 0 0; % negro
    1 1 1;   % blanco
    0 1 0;   % verde
    1 1 0;   % amarillo
    1 0.5 0; % naranja  
    1 0 0;   % rojo
];

%Numero de iteraciones 
it=2500;

contadorAgentes=zeros(1,it+1); %variables para contar cuanto salen del sistema
contadorAgentes(1,1)=NA; %% queda pendiente la explicacion

%3. Se inician las iteraciones y se escriben las reglas
%A.status 0 indica espera
%A.status 1 indica ataque
%A.status 2 indica regresando

for itera=1:it

    if(mod(itera,25) == 0)
        row = randi([12,84]);
        col = randi([47,54]);
        A.meteoritos(row,col) = 1000;
    end

    for h=1:length(M)
        %almacenamiento de las coordenadas originales del agente al inicio de cada iteración.
        posm=M(h);
        posn=N(h);

        Block=A.meteoritos(M(h)-10:M(h)+10,N(h)-10:N(h)+10);
        Lo=ismember(Block,1000);

        % Valida que el agente no tiene un meteorito en su vecindario
        if (sum(Lo(:))>0)
            % si hay un meteorito, la nave pasa a atacar
            A.status(M(h),N(h))=1;
        end

        % si no hay salida, Se evalua el siguiente paso a dar en el mapa
        if (A.status(posm,posn)==2)
            BlockEstacion=D(posm-1:posm+1,posn-1:posn+1);
            minblock = min(min(BlockEstacion));
            if minblock<D(posm,posn)
                [menory,menorx] = find(BlockEstacion == minblock);
                %Se deja el espacio anterior vacio
                A.pos(posm,posn)=0;
                A.status(posm,posn)=0;
                
                if (menorx > 2)
                    posn = posn + 1;
                elseif (menorx < 2)
                    posn = posn - 1;
                end

                if (menory > 2)
                    posm = posm + 1;
                elseif (menory < 2)
                    posm = posm - 1;
                end
                A.status(posm,posn)=2;
                A.pos(posm,posn)=1;
            else
                A.status(posm,posn)=0;
            end
        end
                
        % Si esta en ataque
        if (A.status(posm,posn)==1)
            buscarmeteoritos = find(A.meteoritos == 1000);
            if(buscarmeteoritos)
                [meteorito_pos_y,meteorito_pos_x] = find(A.meteoritos == 1000);
                index = randi(1,length(meteorito_pos_x));
                meteorito_pos_x = meteorito_pos_x(index);
                meteorito_pos_y = meteorito_pos_y(index);

                A.pos(posm,posn)=0;
                A.status(posm,posn)=0;

                if (meteorito_pos_x > posn)
                    posn = posn + 1;
                elseif (meteorito_pos_x < posn)
                    posn = posn - 1;
                end

                if (meteorito_pos_y > posm)
                    posm = posm + 1;
                elseif (meteorito_pos_y < posm)
                    posm = posm - 1;
                end

                if(meteorito_pos_x(1,1) == posn(1,1)) && (meteorito_pos_y(1,1) == posm(1,1))
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
    %meteoritos
    
    % Actualización de Posiciones y Visualización
    T1=M(M>0);
    T2=N(N>0);
       
    M=T1;
    N=T2;
       
    u=length(M);
    dat=randperm(u);
       
    M=M(dat);
    N=N(dat);

    % Etiquetar componentes conectados
    cc = bwconncomp(A.pos);
    
    % Crear la matriz de calor basada en los clusters
    heat_matrix = zeros(num_rows, num_cols);
    
    % Asignar un nivel de calor a cada cluster basado en su tamaño
    for i = 1:cc.NumObjects
        cluster_size = numel(cc.PixelIdxList{i});
        heat_level = min(cluster_size, 5); % Máximo nivel de calor
        heat_matrix(cc.PixelIdxList{i}) = heat_level;
    end

    % Crear una matriz de colores basada en la matriz de calor
    color_matrix = zeros(num_rows, num_cols, 3);
    for i = 1:num_rows
        for j = 1:num_cols
            color_idx = heat_matrix(i, j);
            if color_idx > 0
                color_matrix(i, j, :) = cmap(color_idx + 1, :);
            end
        end
    end
    
    meteoritos = find(A.meteoritos == 1000);
    if(meteoritos)
        for i=1:length(meteoritos)
            [meteorito_y, meteorito_x] = ind2sub(size(D), meteoritos(i)); 
            color_matrix(meteorito_y, meteorito_x, :) = [1, 0, 0]; % Asigna rojo al meteorito
        end
    end
    
    %Se prepara el mapeo del ambiente
    labG =zeros(100,100);
    labG(labBi==255)=1;
    
    %Se grafican los agentes y los muros
    RGB(:,:,3)=labG;

    combined_matrix = max(RGB, color_matrix);
    J = imresize(combined_matrix,5);
    
    figure(1)
    imshow(J);
    hold on
    drawnow;
end

figure(2)
plot(contadorAgentes(1,1:2500),'-','LineWidth', 2)
xlabel('Iterations');
ylabel('Number of agents');
title('Progress of agent evacuation with respect to iterations');
grid on; % Añadir una cuadrícula a la gráfica

toc