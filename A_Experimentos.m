function A_Experimentos()
problems=5;

Ns       = [15,6,19,14,14,4,22,23];
bound_p1 = [0 60;0 60;0 60;0 60;-60 60;-60 60;0 2*pi;-60 60;-60 60;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi];
bound_p2 = [0 50;0 50;0 50;0 50;-50 50;-50 50];
bound_p3 = [0 60;0 60;0 60;0 60;-60 60;-60 60;0 2*pi;-60 60;-60 60;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi;0 2*pi];
bound_p4 = [0 100;0 100;0 100;-50 0;-50 0;-50 0;0 150;0 150;0 50;0 150;0 150;0 150;(pi/12) pi;((3*pi)/4) ((5*pi)/4)];
bound_p5 = [0 100;0 100;0 100;-50 0;-50 0;-50 0;0 150;0 150;0 50;0 150;0 150;0 150;(pi/12) pi;((3*pi)/4) ((5*pi)/4)];

bounds = {bound_p1,bound_p2,bound_p3,bound_p4,bound_p5};
E = [400000,20000, 225000,125000,200000,12000];


[stat,struc] = fileattrib;
PathCurrent = struc.Name;
GlobalFolderName='/Resultados de Experimentos';
mkdir(PathCurrent, GlobalFolderName);

FolderName=['/A'];
PathFolder = [PathCurrent GlobalFolderName FolderName];
mkdir([PathCurrent GlobalFolderName], FolderName);

for i=1:problems
    Problema=i
    Evaluaciones=E(i)
    set_problem_i(i);
    set_problem_n(Ns(i));
    set_problem_bound(bounds{i});
    
    Solutions=[];
    for j=1:3
        corrida=j
        [NM1]= NM(E(i));
        [NM2]= NMS(E(i));
        [NM3]= NMELA(E(i));
        [NM4]= NMEILA(E(i));
        [NM5]= NM2ELA(E(i));
        
        Solutions(j,:)=[NM1(Ns(i)+1),NM2(Ns(i)+1),NM3(Ns(i)+1),NM4(Ns(i)+1),NM5(Ns(i)+1)]
    end
    mejores=min(Solutions);
    peores=max(Solutions);
    medianas=median(Solutions);
    promedios=mean(Solutions);
    desvaciones_estandar=std(Solutions);
    
    for a=1:5
        display(['ALGORTIMO ',num2str(a)]);
        mejor=mejores(a)
        peor=peores(a)
        mediana=medianas(a)
        promedio=promedios(a)
        desviacion=desvaciones_estandar(a)
    end
        % Salvando resultados de Ejecuciones
    Name = ['/A-Resultados_Estadistica_Problema_' num2str(i) '.xls'];
    Nameexcel = [PathFolder Name];
    encabezado = {' ' 'NM' 'NMS'  'NMELA'  'NMEILA'  'NM2ELA'};
    xlswrite(Nameexcel, encabezado, 1, 'B1:F1');
    xlswrite(Nameexcel, Solutions, 1, 'B2:F32');
    xlswrite(Nameexcel, {'Mejor'}, 1, 'A33');
    xlswrite(Nameexcel, mejores, 1, 'B33:F33');
    xlswrite(Nameexcel, {'Peor'}, 1, 'A34');
    xlswrite(Nameexcel, peores, 1, 'B34:F34');
    xlswrite(Nameexcel, {'Mediana'}, 1, 'A35');
    xlswrite(Nameexcel, medianas, 1, 'B35:F35');
    xlswrite(Nameexcel, {'Promedio'}, 1, 'A36');
    xlswrite(Nameexcel, promedios, 1, 'B36:F36');
    xlswrite(Nameexcel, {'Desv. Stan'}, 1, 'A37');
    xlswrite(Nameexcel, desvaciones_estandar, 1, 'B37:F37');
 
    % Pruebas de Firedman y Bonferroni 
    disp('******* Friedman Test results ********')
    [p,tbl,stats] = friedman(Solutions,1)
    disp('******* Bonferrony Test results ******')
    [results,means,h]= multcompare(stats,'CType','bonferroni')
    % Salvando resultados de prueba de Friedman
    Name = ['/A-Resultados_Firedman_Problema_' num2str(i) '.xls'];
    Nameexcel = [PathFolder Name];
    xlswrite(Nameexcel, encabezado, 1, 'B1:F1');
    xlswrite(Nameexcel, {'Ranks'}, 1, 'A2');
    xlswrite(Nameexcel, stats.meanranks, 1, 'B2:F2');    
    xlswrite(Nameexcel, {'P-value'}, 1, 'A3');
    xlswrite(Nameexcel, p, 1, 'B3');
    % Salvando gráficas de Bonferroni 
    saveas(h,fullfile(PathFolder, ['Bonferroni_HNMED_VS_ED' num2str(i)]), 'fig');
    saveas(h,fullfile(PathFolder, ['Bonferroni_HNMED_VS_ED' num2str(i)]), 'png');
    beep;
end

