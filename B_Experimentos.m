function B_C_Experimentos()
clear
problems=6;

Ns       = [15,6,19,14,14,4];
bound_p1 = [0 60;0 60;0 60;0 60;-60 60;-60 60;0 2*pi;-60 60;-60 60;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi];
bound_p2 = [0 50;0 50;0 50;0 50;-50 50;-50 50];
bound_p3 = [0 60;0 60;0 60;0 60;-60 60;-60 60;0 2*pi;-60 60;-60 60;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi];
bound_p4 = [0 100;0 100;0 100;-50 0;-50 0;-50 0;0 150;0 150;0 50;0 150;0 150;0 150;(pi/12) pi;((3*pi)/4) ((5*pi)/4)];
bound_p5 = [0 100;0 100;0 100;-50 0;-50 0;-50 0;0 150;0 150;0 50;0 150;0 150;0 150;(pi/12) pi;((3*pi)/4) ((5*pi)/4)];

bounds = {bound_p1,bound_p2,bound_p3,bound_p4,bound_p5};
E = [750030,5000, 450018,750030,450018,12000];
S = [9,  5,   7,   9,   9   ,7];
NP = [138,   56,   138,   138,   138   ,50];


[stat,struc] = fileattrib;
PathCurrent = struc.Name;
GlobalFolderName='/Resultados de Experimentos';
mkdir(PathCurrent, GlobalFolderName);

FolderName=['/B'];
PathFolder = [PathCurrent GlobalFolderName FolderName];
mkdir([PathCurrent GlobalFolderName], FolderName);

FolderName2=['/B'];
PathFolder2 = [PathCurrent GlobalFolderName FolderName2];
mkdir([PathCurrent GlobalFolderName], FolderName2);

for i=1:problems
    Problema=i
    Evaluaciones=E(i)
    N=Ns(i)
    Cantidad_Simplex=S(i)
    Poblacion=S(i)*(Ns(i)+1)
    set_problem_i(i);
    set_problem_n(Ns(i));
    set_problem_bound(bounds{i});
    
    Solutions=[];
    for j=1:3
        corrida=j
        if i~=6
        [DErand1Bin, E_F1]=DE_rand_1_bin(E(i),NP(i));
        [HNMED1,E_F2]= Nelder_Mead_ED1(E(i),S(i));
        [HNMED2,E_F3]= Nelder_Mead_ED2(E(i),S(i));
        [HNMED3,E_F4]= Nelder_Mead_ED3(E(i),S(i));
        [HNMED4,E_F5]= Nelder_Mead_ED4(E(i),S(i));
        [HNMED5,E_F6]= Nelder_Mead_ED5(E(i),S(i));
        else
        %___________________Micro-Red______________________________________
        % El problema de la microred electrica se implemta diferente llamar
        % aqui los métodos para el problema
        [DErand1Bin, E_F1]=MR_DE_rand_1_bin(E(i),50,j);
        [HNMED1, E_F2]= MR_HNMED1(E(i),S(i),j);
        [HNMED2, E_F3]= MR_HNMED2(E(i),S(i),j);
        [HNMED3, E_F4]= MR_HNMED3(E(i),S(i),j);
        [HNMED4, E_F5]= MR_HNMED4(E(i),S(i),j);
        [HNMED5, E_F6]= MR_HNMED5(E(i),S(i),j);
        
        end
        Solutions(j,:)=[DErand1Bin(Ns(i)+1),HNMED1(Ns(i)+1),HNMED2(Ns(i)+1),HNMED3(Ns(i)+1),HNMED4(Ns(i)+1),HNMED5(Ns(i)+1)]
        % Arreglos con tasas de convergencia para cada ejecución
        E_F1s{j}=E_F1(1:end,:,:);
        E_F2s{j}=E_F2(1:end,:,:);
        E_F3s{j}=E_F3(1:end,:,:);
        E_F4s{j}=E_F4(1:end,:,:);
        E_F5s{j}=E_F5(1:end,:,:);
        E_F6s{j}=E_F6(1:end,:,:);
    end
    
    % Hallar los índices correspondientes a la media de las muestras
    media_index=[];
    s=size(Solutions);
    for l=1:s(2)
        C=Solutions(:,l);
        index=find(C==median(C));
        if(length(index)>0)
            media_index(l)=index(1);
        else
            media_index(l)= index;
        end
    end
    
    % Seleccionar la tasa de convergencia correspondiente a la media de
    % cada muestra
    E_F1=E_F1s{media_index(1)};
    E_F2=E_F2s{media_index(2)};
    E_F3=E_F3s{media_index(3)};
    E_F4=E_F4s{media_index(4)};
    E_F5=E_F5s{media_index(5)};
    E_F6=E_F6s{media_index(6)};
    % Grafica funciones de convergencia.
    figure();
    q=plot(E_F1(:,2),E_F1(:,1),E_F2(:,2),E_F2(:,1),E_F3(:,2),E_F3(:,1),E_F4(:,2),E_F4(:,1),E_F5(:,2),E_F5(:,1),E_F6(:,2),E_F6(:,1));
    legend('ED','HNMED-V1','HNMED-V2','HNMED-V3','HNMED-V4','HNMED-V5')
    set(q,'linewidth',1);
    % Para el caso de la micro-red usar escala lineal para valores de la 
    % FO negativos 
    if i~=6
    set(gca,'YScale','log');
    set(gca,'XScale','log');
    else
    set(gca,'YScale','linear');
    set(gca,'XScale','linear');
    end
    title('Gráfica f(x) vs NEF')
    xlabel('NEF');
    ylabel('f(x)');
    % Salvar gráfica
    saveas(gcf,fullfile(PathFolder, ['Grafica_Convergencia_Problema_' num2str(i)]), 'png');
    saveas(gcf,fullfile(PathFolder, ['Grafica_Convergencia_Problema_' num2str(i)]), 'fig');
    mejores=min(Solutions);
    peores=max(Solutions);
    medianas=median(Solutions);
    promedios=mean(Solutions);
    desvaciones_estandar=std(Solutions);
    
    for a=1:6
        display(['ALGORTIMO ',num2str(a)]);
        mejor=mejores(a)
        peor=peores(a)
        mediana=medianas(a)
        promedio=promedios(a)
        desviacion=desvaciones_estandar(a)
    end
        % Salvando resultados de Ejecuciones
    Name = ['/Resultados_Estadistica_Problema_' num2str(i) '.xls'];
    Nameexcel = [PathFolder Name];
    encabezado = {' ' 'ED' 'HNMED1'  'HNMED2'  'HNMED3'  'HNMED4'  'HNMED5'};
    xlswrite(Nameexcel, encabezado, 1, 'B1:G1');
    xlswrite(Nameexcel, Solutions, 1, 'B2:G32');
    xlswrite(Nameexcel, {'Mejor'}, 1, 'A33');
    xlswrite(Nameexcel, mejores, 1, 'B33:G33');
    xlswrite(Nameexcel, {'Peor'}, 1, 'A34');
    xlswrite(Nameexcel, peores, 1, 'B34:G34');
    xlswrite(Nameexcel, {'Mediana'}, 1, 'A35');
    xlswrite(Nameexcel, medianas, 1, 'B35:G35');
    xlswrite(Nameexcel, {'Promedio'}, 1, 'A36');
    xlswrite(Nameexcel, promedios, 1, 'B36:G36');
    xlswrite(Nameexcel, {'Desv. Stan'}, 1, 'A37');
    xlswrite(Nameexcel, desvaciones_estandar, 1, 'B37:G37');
 
    % Pruebas de Firedman y Bonferroni 
    disp('******* Friedman Test results ********')
    [p,tbl,stats] = friedman(Solutions,1)
    disp('******* Bonferrony Test results ******')
    [results,means,h]= multcompare(stats,'CType','bonferroni')
    % Salvando resultados de prueba de Friedman
    Name = ['/Resultados_Firedman_Problema_' num2str(i) '.xls'];
    Nameexcel = [PathFolder2 Name];
    encabezado = {' ' 'ED' 'HNMED1'  'HNMED2'  'HNMED3'  'HNMED4'  'HNMED5'};
    xlswrite(Nameexcel, encabezado, 1, 'B1:G1');
    xlswrite(Nameexcel, {'Ranks'}, 1, 'A2');
    xlswrite(Nameexcel, stats.meanranks, 1, 'B2:G2');    
    xlswrite(Nameexcel, {'P-value'}, 1, 'A3');
    xlswrite(Nameexcel, p, 1, 'B3');
    % Salvando gráficas de Bonferroni 
    saveas(h,fullfile(PathFolder2, ['Bonferroni_HNMED_VS_ED' num2str(i)]), 'fig');
    saveas(h,fullfile(PathFolder2, ['Bonferroni_HNMED_VS_ED' num2str(i)]), 'png');
    beep;
end

end