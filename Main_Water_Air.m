% Programa principal parte 1 para tratamento de imagem no Rotor BCS 
% Imagens providas pela camera nova 10/10/2018
% Marcos Mendes - 144124
% ----------- AGUA E AR -----------
clc; close all; clearvars;
%retirando o warning das imagens
warning('off', 'Images:initSize:adjustingMag'); 

% abertura do video avi
im=3; %Para mudar o arquivo de abertura ,,
lente=50; %quando a lente = 50 existe um filtro morfologico carregado para remover as sombras
if lente == 50
    eccentr = .96; % a lente de 50 mm apresenta mais sombras, que envolve a bolha. O filtro 
    % peculiar precisa ser mais "rigido" nesse aspecto (AJUSTAR) 
else
    eccentr=.96;
end

% nome=['C:\Users\user-pc\Videos\900rpm_0012g\im',num2str(im),'_',num2str(lente),'mm.avi']; 
nome=['C:\Users\user-pc\Videos\900rpm_0025g\im',num2str(im),'_',num2str(lente),'mm.avi']; 
Bifasico = VideoReader(nome);

clear nome 
% Dados do vídeo:altura, comprimento e número de frames
I = Bifasico.Height;
J = Bifasico.Width;

I1 = read(Bifasico,1);
figure(1)
imshow(I1)

% gatilhos do programa
maskok=0;
Omega= 900; %rotacao da bomba (600, 900, 1200)
Fr=1500; %frequencia de aquisicao (1000, 1500, 2000)
opcao=1; %corte da mascara: = 1 para deixar apenas o canal =2 para deixar o impelidor todo
filtro=1; % =1 para ver as imagens criadas dentro da funcao filtro; =2 para nao ver
canal=1; % define a quantidade de canais que sera usada no tratamento de imagem

% 1 - criacao da regiao de interesse - MASCARA
if maskok==1
    file= ['Masks\',num2str(Omega),'\im',num2str(im),'.mat'];
    mask=zeros(J,I,'double');%alocação do frame Gray uint8;
    for n = 1:canal;
        figure(1);clf;imshow(I1)
        [mask1, xi, yi] = roipoly(I1);
        %atualiza a máscara
        mask = mask+mask1;
    end
    save(file, 'mask', 'xi', 'yi')
else
    file= ['Masks\',num2str(Omega),'\im',num2str(im),'.mat'];
    load(file);
end

mask1=mask;
[x,y,~,~]=checkmask(mask1,1); %funcao para checar a posicao da mascara; x2 e y2 mascara toda
clear canal maskok mask x2 y2 
% Inicializacao das matrizes cubicas
centroids_cb=NaN.*zeros(50,2,205); % centroide
diamAreas_cb=NaN.*zeros(70,2,205); % diametros
myCell=cell(4,205);
size_cb=NaN.*zeros(1,1,205); % tamanho de linhas da matriz centroide
thre= 0.37; %This threshold is related to the image difference
% im1= 0.40 ; im3 = 0.41  ; im5= 0.40  ; im7 = 0.41  ; mg0006 - 900rpm
% im1= 0.42 ; im3 = 0.47  ; im5= 0.43  ; im7 = 0.44  ; mg0012 - 900rpm
% im1= 0.36 ; im3 = 0.37  ; im5= 0.43  ; im7 = 0.43  ; mg0025 - 900rpm
% im1= 0.44 ; im3 = 0.42  ; im5= 0.42  ; im7 = 0.42  ; mg0025 - 600rpm
% im1= 0.37 ; im3 = 0.38  ; im5= 0.37  ; im7 =   ; mg0025 - 1200rpm

thre2= 0.31; %This threshold is related to imaget1 and imaget2
% im1= 0.22 ; im3 = 0.24  ; im5= 0.24  ; im7 = 0.21  ; mg0006 - 900rpm
% im1= 0.35 ; im3 = 0.36  ; im5= 0.36  ; im7 = 0.29  ; mg0012 - 900rpm
% im1= 0.16 ; im3 = 0.31  ; im5= 0.18  ; im7 = 0.18  ; mg0025 - 900rpm
% im1= 0.18 ; im3 = 0.28  ; im5= 0.14  ; im7 = 0.18  ; mg0025 - 600rpm
% im1= 0.17 ; im3 = 0.14  ; im5= 0.22  ; im7 =   ; mg0025 - 1200rpm

for i = 1 : 250
    
    I1 = read(Bifasico, i);
    I2 = read(Bifasico, i+1);
    
    % ROTACIONANDO A IMAGEM (rever apos aplicar os filtros) 
    J1=imrotate(I1,[(i-1)*(Omega/60*(360/Fr))],'bilinear','crop');
    J2=imrotate(I2,[i*(Omega/60*(360/Fr))],'bilinear','crop');
    
    % CORTANDO A IMAGEM
    % Recortando com base na mascara
    CI_withmask1 = J1; CI_withmask2 = J2; %if the image was in RGB, convert here. 
    clear J1 J2 J3
    CI_withmask1(~mask1) = 0; CI_withmask2(~mask1) = 0; % It is not necessary rotate the mask because J1 is already been rotated
    
    % disp('Recorte da região de interesse');
    [CI1]=crop(I1,1); [CI2]=crop(I2,1);%Variaveis: (imagem; regiao onde sera cortada 1 == externa, 2 == centro)
    [CI1]=crop(CI1,2); [CI2]=crop(CI2,2); 
    
    % Recortando apenas 1 canal, de modo retangular
    [CI1,mask]=crop_canal(CI_withmask1,x,y,opcao,mask1); [CI2,~]=crop_canal(CI_withmask2,x,y,opcao,mask1);
    [x2,y2,~,~]=checkmask(mask,2); %Obtendo a transformada da mascara para reverter depois no sistema do rotor  save xy_mg0006.mat x2 y2
    
    if i == 1
        figure(2) %mostrando a regiao de trabalho para o operador
%         subplot(1,2,1)
        imshow(CI1)
%         subplot(1,2,2)
%         imshow(mask)
    end
    % FILTROS
    [CI1f4]=filtros2(CI1,mask,2,4); % (~,~,~,k) o valor k vai de 1 a 7
%     fig(:,:,i)=CI1f4; 
    
    % Morph Operation => Removing shadows
    if lente == 50
       [CI1]=filtros2(CI1,lente,2,11); 
       [CI2]=filtros2(CI2,lente,2,11); 
    end
    
    %% IMAGE SUBTRACTION
    [CI1f4]=filtros2(CI1,mask,2,4);
    [CI2f4]=filtros2(CI2,mask,2,4);

    Kf4 = imabsdiff(CI1f4,CI2f4);
    [Kf4f4]=filtros2(Kf4,mask,2,4);
    Kf4f4mf=medfilt2(Kf4f4,[3 3],'symmetric');
    
    if i == 1
        figure(3), imshow(Kf4)
        title('Kf4', 'FontSize', 16);
        figure(4), imshow(Kf4f4)
        title('Kf4f4 - with adjust filter applied', 'FontSize', 16);
        figure(5), imshow(Kf4f4mf)
        title('Kf4f4mf - with median filter applied', 'FontSize', 16);       
    end
    
    [CI1f5]=filtros2(CI1f4,mask,2,10); [CI1f5]=filtros2(CI1f5,mask,2,5);
    [CI2f5]=filtros2(CI2f4,mask,2,10); [CI2f5]=filtros2(CI2f5,mask,2,5);
    CI1f5(~mask) = 0; CI2f5(~mask) = 0;
   
%     %Nao ficou bom
%     % -------------
%     Kf5 = imabsdiff(CI1f5,CI2f5);
%     figure, imshow(Kf5)    
%     
%     %Nao ficou bom
%     % -------------
%     [CI1f6]=filtros2(CI1f4,mask,2,6);
%     [CI2f6]=filtros2(CI2f4,mask,2,6);
%     
%     Kf6 = imabsdiff(CI1f6,CI2f6);
%     figure, imshow(Kf6)
%     title('Kf6', 'FontSize', 16);
%     
%     %Nao ficou bom
%     % -------------
%     [CI1f8]=filtros2(CI1,mask,2,8);
%     [CI2f8]=filtros2(CI2,mask,2,8);
%     
%     Kf8 = imabsdiff(CI1f8,CI2f8);
%     figure, imshow(Kf8)
%     title('Kf8', 'FontSize', 16);
     
    
%%    BINARIZATION
%     Inicializacao do tratamento da imagem com filtros
%     ----------------------------------------
%     usando o ajuste da imagem e comparando

    K_adjust=imadjust(Kf4f4); 
    [K_adjust]=filtros2(K_adjust,mask,2,10); %atencao para o ajuste de contraste
    
    levelK=graythresh(Kf4f4);
    levelK_adjust=graythresh(K_adjust);
    clear levelCI_adjust1 levelCI_adjust2 
    
    if i == 1
        figure(15)
        subplot(2,2,1)
        imshow(Kf4f4)
        title(levelK, 'FontSize', 16);
        subplot(2,2,2)
        imhist(Kf4f4)
        subplot(2,2,3)
        imshow(K_adjust)
        title(levelK_adjust, 'FontSize', 16);
        subplot(2,2,4)
        imhist(K_adjust)
    end   
    
    % Average, median filter and wiener to remove noise
    Kaverage = filter2(fspecial('average',3),K_adjust)/255;
    Kmedian = medfilt2(K_adjust,[4 4],'symmetric');
    [Kwiener,kwienoise] = wiener2(K_adjust,[4 4]);
    KmedianWiener = medfilt2(Kwiener,[5 5],'symmetric');

    if i == 1
        figure(7), imshowpair(Kaverage,Kmedian,'montage')
        title('Kaverage and Kmedian', 'FontSize', 16);        
        
        figure(8), imshowpair(Kwiener,KmedianWiener,'montage')
        title('KWiener and KmedianWiener', 'FontSize', 16);
    end
    
%     BINARIZACAO DA IMAGEM (binary function)
%     definindo funcoes para binarizacao   
%     f1=@(img) imbinarize(img,'adaptive','ForegroundPolarity','dark','Sensitivity',thre); 
    f1=@(img) imbinarize(img,(levelK+thre)); %0.58 se usar sem o filtro case 10 / 0.47 para um canal
 
    bWadj = roifilt2(K_adjust,mask,f1); bWadj(~mask) = 0;
    bWadj = medfilt2(bWadj,[5 5],'symmetric');
    
    bWave = roifilt2(Kaverage,mask,f1); bWave(~mask) = 0;
    bWave = medfilt2(bWave,[5 5],'symmetric');
    
    bWmed = roifilt2(Kmedian ,mask,f1); bWmed(~mask) = 0; 
    bWmed=bwareaopen(bWmed, 5);
    
    bWwie = roifilt2(Kwiener,mask,f1); bWwie(~mask) = 0;
    bWwie = medfilt2(bWwie,[5 5],'symmetric'); 
    
    bWmed1 = Kmedian < 160; %Imagem binaria; 
    bWmed1 = bwareaopen(~bWmed1, 10);  % Esse parte é uma opcao para remover area menor que 10
    bWmed2 = Kmedian < 190; %Imagem binaria; 
    bWmed3 = Kmedian < 195; %Imagem binaria; 
    bWmed4 = Kmedian < 200; %Imagem binaria;    

    if i == 1
        figure(6)
        imshow(~bWmed)
        title(sprintf('Median ft, %2f', std(im2double(bWmed(:)))), 'FontSize', 12);
        
        figure(9)
        subplot(2,2,1)
        imshow(~bWadj)
        title(sprintf('Img orig, %2f', std(im2double(bWadj(:)))), 'FontSize', 12);
        subplot(2,2,2)
        imshow(~bWave)
        title(sprintf('Average ft, %2f', std(im2double(bWave(:)))), 'FontSize', 12);
        subplot(2,2,3)
        imshow(~bWmed)
        title(sprintf('Median ft, %2f', std(im2double(bWmed(:)))), 'FontSize', 12);
        subplot(2,2,4)
        imshow(~bWwie)
        title(sprintf('Wiener ft, %2f', std(im2double(bWwie(:)))), 'FontSize', 12);
        
        figure(10)
        subplot(2,2,1)
        imshow(~bWmed1)
        title(sprintf('Img < 160, %2f', std(im2double(bWmed1(:)))), 'FontSize', 12);
        subplot(2,2,2)
        imshow(bWmed2)
        title(sprintf('Img < 190, %2f', std(im2double(bWmed2(:)))), 'FontSize', 12);
        subplot(2,2,3)
        imshow(bWmed3)
        title(sprintf('Img < 195, %2f', std(im2double(bWmed3(:)))), 'FontSize', 12);
        subplot(2,2,4)
        imshow(bWmed4)
        title(sprintf('Img < 200, %2f', std(im2double(bWmed4(:)))), 'FontSize', 12);
    end 
    
    clear bWmed1 bWmed2 bWmed3 bWmed4 bWave bWadj bWwie f1 f2 Kaverage Kmedian KmedianWiener Kwiener Kwienernoise levelK levelK_adjust ...
        signal_var Kf4 Kf4f4
    % Aparentemente a mediana foi a melhor. Escolhi a da figura 17 pois a
    % variação do threshold sera importante, ja que cada canal vai apresentar variacao. 
    bWmed = ~bWmed; bWmed(~mask) = 0;
    % (binary function end)
    
    %% TRATAMENTO FINAL DAS IMAGENS INDIVIDUAIS (C1 e C2) 
    if lente == 50
        [bWCI1med]=binary(CI1f4,i,mask,thre2,lente); 
%         bWCI1med = imfill(~bWCI1med,'holes'); bWCI1med =~bWCI1med;
        [bWCI2med]=binary(CI2f4,i,mask,thre2,lente); 
%         bWCI2med = imfill(~bWCI2med,'holes'); bWCI2med =~bWCI2med;
    else
    % O Laplaciano parece estar funcionando melhor para lente de 80 mm!
        [bWCI1med]=binary(CI1f5,i,mask,thre2,lente); 
        [bWCI2med]=binary(CI2f5,i,mask,thre2,lente); 
    end
    %% STEP 1 AND 2
    % 1: Get point to be excluded in next step 
    % ----------------------------------------------------------------
    % O filtro "peculiar" ira encontrar peculiaridades, tal como area e
    % excentricidade fora do esperado
    [L2,B,A,uK]=peculiar(bWmed,eccentr);
    [centroidsK, ~, ~, ~] = variablesWA(bWmed,uK);
    % O filtro "variablesWA" ira encontrar os centroides, diametros e areas
    
    % 2: Excluded some lines in the matrix that are out from specifications
    % ----------------------------------------------------------------
    [l1,b1,a1,u1]=peculiar(bWCI1med,eccentr);
    [centroids1, diameters1, area1, ~] = variablesWA(bWCI1med,u1);
    [~,~,~,u2]=peculiar(bWCI2med,eccentr);
    [centroids2, diameters2, area2, white] = variablesWA(bWCI2med,u2); 
 
    if  i < 10
        Data_test(L2,A,B,uK,18) %verificando o que sera removido da diferenca das imagens
        Data_test(l1,a1,b1,u1,19) %verificando o que sera removido das imagens
    end
    
    clear L2 B A l1 b1 a1 
    %% STEP 3: Comparing the displacement matrix (K) with 2 timesteps from move matrix
    % ----------------------------------------------------------------
    % Contando os elementos para saber se a quantidade de bolhas detectada
    % é a mesma nos dois passos de tempo considerados 
    [kK,~]=size(centroidsK); % contando a diferenca
    [k1,~]=size(centroids1); % contando o frame atual
    [k2,~]=size(centroids2); % contando o frame posterior    
    
    [auxiliar1,auxiliarK,diameters1,area1] = pareaStep1WA(centroids1,centroidsK,16,16,k1,kK,diameters1,area1); %(A e a menor matriz
    [auxiliar2,~,diameters2,area2] = pareaStep1WA(centroids2,centroidsK,16,16,k2,kK,diameters2,area2);
    clear kK
    if i < 5
        figure(21)
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
 
        subplot(1,2,1)
        scatter(centroidsK(:,1),centroidsK(:,2),'g+')
        hold on
        scatter(centroids1(:,1),centroids1(:,2),'ro')
        scatter(centroids2(:,1),centroids2(:,2),'bs')
        plot(x2,y2,'k')
        set(gca,'Ydir','reverse')
%         title('Without Comparison','FontSize', 16);
        axis equal
        ylim([-10 480]); 
        xlim([0 750]);
        set(gca,'FontSize',16);
        
        subplot(1,2,2)
        scatter(centroidsK(:,1),centroidsK(:,2),'g+')
        hold on
        scatter(auxiliar1(:,1),auxiliar1(:,2),'ro')
        scatter(auxiliar2(:,1),auxiliar2(:,2),'bs')
        plot(x2,y2,'k')
        set(gca,'Ydir','reverse')
%         title('After Comparison','FontSize', 16);
        axis equal
        ylim([-10 480]); 
        xlim([0 750]);
        legend('Subtraction','First frame', 'Second frame', 'Location', 'southwest' );
        set(gca,'FontSize',16);
    end
    [k1,~]=size(auxiliar1);
    
    if i == 1
        figure(22) %Figura que compara as imagens da bolha em t1 e t2
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        subplot(1,2,1)
        imshow(CI1f4); 
%         scatter(auxiliar1(:,1),auxiliar1(:,2),'ro')
        title('CI1f4', 'FontSize', 16);
        subplot(1,2,2)
        imshow(CI2f4); 
%         scatter(auxiliar2(:,1),auxiliar2(:,2),'bs')
        title('CI2f4', 'FontSize', 16);
    end
    if i == 1
        figure(23);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        subplot(1,2,1)
        imshow(bWCI1med); 
        title('CI1f4', 'FontSize', 16);
        subplot(1,2,2)
        imshow(bWCI2med); 
        title('CI2f4', 'FontSize', 16);
    end    
    %% STEP 4: Saving on cubic matrix
    % ---------------------------------------------------------------- 
    % correcting the system change (because the code operates in one
    % channel) figure, imshow(CI_withmask1)
    auxiliar1(:,1)=auxiliar1(:,1)+(min(x)-1);auxiliar1(:,2)=auxiliar1(:,2)+(min(y)-1);
    if i<3 %verificando a rotacao/mascara
        figure(24);
        imshow(CI_withmask1)
        hold on
        scatter(auxiliar1(:,1),auxiliar1(:,2),'ro'); hold off       
    end
    
    % saving
    centroids_cb(1:k1,:,i)=auxiliar1(1:k1,:); % centroide
    diamAreas_cb(1:k1,1,i)=diameters1; % diametros 
    diamAreas_cb(1:length(area1),2,i)= area1;% areas
    size_cb(1,1,i)=k1; % 
    myCell(:,i) = {auxiliar1(1:k1,:); diameters1; area1;k1};    
        
end

disp('finish') %aqui vc deve salvar save mg0025_im.mat centroids_cb diamAreas_cb size_cb
file2= ['Masks\',num2str(Omega),'\im',num2str(im),'Result.mat'];
save(file2, 'myCell')
% save mg0025_im.mat myCell
