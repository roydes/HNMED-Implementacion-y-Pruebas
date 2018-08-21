function F_Experimentos()
clear
problems=6;

Ns       = [15,6,19,14,14,4];
bound_p1 = [0 60;0 60;0 60;0 60;-60 60;-60 60;0 2*pi;-60 60;-60 60;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi];
bound_p2 = [0 50;0 50;0 50;0 50;-50 50;-50 50];
bound_p3 = [0 60;0 60;0 60;0 60;-60 60;-60 60;0 2*pi;-60 60;-60 60;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi];
bound_p4 = [0 100;0 100;0 100;-50 0;-50 0;-50 0;0 150;0 150;0 50;0 150;0 150;0 150;(pi/12) pi;((3*pi)/4) ((5*pi)/4)];
bound_p5 = [0 100;0 100;0 100;-50 0;-50 0;-50 0;0 150;0 150;0 50;0 150;0 150;0 150;(pi/12) pi;((3*pi)/4) ((5*pi)/4)];

bounds   = {bound_p1,bound_p2,bound_p3,bound_p4,bound_p5};
E_SHADE  = [400000,15000,200000,325000,150000,10000];
E_HNMED5 = [400000,15000, 250000,125000,250000,9000];

S_HNMED5 = [7,  4,   7,   9,   9   ,7];
S_HNMED6 = [9,  4,   7,   9,   9   ,7];
NP = [138,   50,   138,   138,   138   ,50];


[stat,struc] = fileattrib;
PathCurrent = struc.Name;
GlobalFolderName='/Resultados de Experimentos';
mkdir(PathCurrent, GlobalFolderName);

FolderName=['/D'];
PathFolder = [PathCurrent GlobalFolderName FolderName];
mkdir([PathCurrent GlobalFolderName], FolderName);

FolderName2=['/D'];
PathFolder2 = [PathCurrent GlobalFolderName FolderName2];
mkdir([PathCurrent GlobalFolderName], FolderName2);
sample_size=31
for i=1:problems
    Problema=i
    Evaluaciones=E_HNMED6(i)
    N=Ns(i)
    Cantidad_Simplex=S_HNMED6(i)
    Poblacion=S_HNMED6(i)*(Ns(i)+1)
    set_problem_i(i);
    set_problem_n(Ns(i));
    if i~=6
    set_problem_bound(bounds{i});
    end
    Solutions=[];
    Vectors=[];
    for j=1:sample_size
       
        corrida=j
        if i~=6
        [C_LSHADE,E_F1]=CLSHADE(Ns(i),bounds{i},E_SHADE(i),j);    
        [HNMED5,E_F2]= Nelder_Mead_ED5(E_HNMED5(i),S_HNMED5(i));
        C_LSHADE
        HNMED5
        Vectors(j,1,1:Ns(i)+2)=C_LSHADE;
        Vectors(j,2,1:Ns(i)+2)=HNMED5;
        Solutions(j,:)=[C_LSHADE(Ns(i)+1),HNMED5(Ns(i)+1)]

        else
        %___________________Micro-Red______________________________________
        % El problema de la microred electrica se implemta diferente llamar
        % aqui los métodos para el problema
        [C_LSHADE,E_F1]= MR_CLSHADE(j);
        [HNMED5,E_F2]= MR_HNMED5 (E_HNMED5(i),S_HNMED5(i),j);
        Solutions(j,:)=[C_LSHADE, HNMED6]

        end
        % Arreglos con tasas de convergencia para cada ejecución
        E_F1s{j}=E_F1(1:end,:,:);
        E_F2s{j}=E_F2(1:end,:,:);
    end
    mejores=min(Solutions);
    peores=max(Solutions);
    medianas=median(Solutions);
    promedios=mean(Solutions);
    desvaciones_estandar=std(Solutions);
    
    for a=1:2
        display(['ALGORTIMO ',num2str(a)]);
        mejor=mejores(a)
        peor=peores(a)
        mediana=medianas(a)
        promedio=promedios(a)
        desviacion=desvaciones_estandar(a)
    end
    
    % Hallar los índices correspondientes a la media de las muestras
    media_index=[];
    min_index=[];
    s=size(Solutions);
    for l=1:s(2)
        C=Solutions(:,l);
        index_med=find(C==median(C));
        index_min=find(C==min(C));
        if(length(index_med)>0)
            media_index(l)=index_med(1);
        else
            media_index(l)= index_med;
        end
        if(length(index_min)>0)
            min_index(l)=index_min(1);
        else
            min_index(l)= index_min;
        end
    end
    
    % Seleccionar la tasa de convergencia correspondiente a la media de
    % cada muestra
    E_F1=E_F1s{media_index(1)};
    E_F2=E_F2s{media_index(2)};
    
    %*************** SALVANDO EXPERMIENTO F   ****************************
    % Grafica funciones de convergencian Experimento.
    figure();
    q=plot(E_F1(:,2),E_F1(:,1),E_F2(:,2),E_F2(:,1));
    legend('C-LSHADE','HNMED-V5')
    set(q,'linewidth',2);
    set(q, {'color'}, {'k'; [51, 51, 204]/255});
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
    saveas(gcf,fullfile(PathFolder, ['F-Grafica_Convergencia_Problema_' num2str(i)]), 'png');
    saveas(gcf,fullfile(PathFolder, ['F-Grafica_Convergencia_Problema_' num2str(i)]), 'fig');
    % Salvando resultados de Ejecuciones
    Name = ['/F-Resultados_Estadistica_Problema_' num2str(i) '.xls'];
    Nameexcel = [PathFolder Name];
    encabezado = {'C-LSHADE'  'HNMED5'};
    xlswrite(Nameexcel, encabezado, 1, 'B1:C1');
    xlswrite(Nameexcel, [Solutions(:,1),Solutions(:,2)], 1, 'B2:C32');
    xlswrite(Nameexcel, {'Mejor'}, 1, 'A33');
    xlswrite(Nameexcel, [mejores(1),mejores(2)], 1, 'B33:C33');
    xlswrite(Nameexcel, {'Peor'}, 1, 'A34');
    xlswrite(Nameexcel, [peores(1),peores(2)], 1, 'B34:C34');
    xlswrite(Nameexcel, {'Mediana'}, 1, 'A35');
    xlswrite(Nameexcel, [medianas(1),medianas(2)], 1, 'B35:C35');
    xlswrite(Nameexcel, {'Promedio'}, 1, 'A36');
    xlswrite(Nameexcel, [promedios(1),promedios(2)], 1, 'B36:C36');
    xlswrite(Nameexcel, {'Desv. Stan'}, 1, 'A37');
    xlswrite(Nameexcel, [desvaciones_estandar(1),desvaciones_estandar(2)], 1, 'B37:C37');
    xlswrite(Nameexcel, {'Mejor indivudio HNMED5'}, 1, 'A38'); 
    best_HNMED5=squeeze(Vectors(min_index(2),2,1:Ns(i)+2));
    xlswrite(Nameexcel,best_HNMED5, 1, ['A39:A' num2str(40+Ns(i)+1)] ); 
    
    
    disp('******* Wilcoxon Test results for Experiment F********')
%   ranksum(x,y) also returns a logical value indicating the test decision.
%   The result h = 1 indicates a rejection of the null hypothesis, and
%   h = 0 indicates a failure to reject the null hypothesis at the 5% 
%   significance level.
%   Right-tailed hypothesis test, where the alternative hypothesis states
%   that the median of x is greater than the median of y.     
     [p,h,stats] =ranksum(Solutions(:,1),Solutions(:,2),'tail','right');
    p_value=p
    h
    rank_sum=stats.ranksum
    winner='HNMED5'
    H0='Aceptada';
    if h
    H0='Rechazada'  ; 
    else
    [p,h,stats] =ranksum(Solutions(:,1),Solutions(:,2),'tail','left');
    if h
        winner='CLSHADE'
    else
       winner='EQUALS' 
    end
    end
    H0
    xlswrite(Nameexcel,{'H0 de Prueba de Wilcoxon'}, 1, 'D1' );
    xlswrite(Nameexcel,{H0}, 1, 'D2' );
    xlswrite(Nameexcel,{'Suma de jerarquías'}, 1, 'E1' );
    xlswrite(Nameexcel,rank_sum, 1, 'E2' );
    xlswrite(Nameexcel,{'Ganador'}, 1, 'E3' );
    xlswrite(Nameexcel,{winner}, 1, 'E4' );
    beep;
end

end