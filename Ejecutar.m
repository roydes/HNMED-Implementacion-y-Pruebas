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
E_HNMED6 = [350000,10000, 175000,125000,125000,9000];
S_HNMED6 = [7,  4,   5,   7,   7   ,2];


[stat,struc] = fileattrib;
PathCurrent = struc.Name;
GlobalFolderName='/Resultados de Experimentos';
mkdir(PathCurrent, GlobalFolderName);

FolderName=['/F'];
PathFolder = [PathCurrent GlobalFolderName FolderName];
mkdir([PathCurrent GlobalFolderName], FolderName);

sample_size=30;
for i=2:2
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
          HJ=EDBOP(E_SHADE(i),35);  
       % [HNMED5,E_F1]= Nelder_Mead_ED5(E_HNMED6(i),S_HNMED6(i));    
        %[C_LSHADE,E_F1]=CLSHADE(Ns(i),bounds{i},E_SHADE(i),j);    
        %[HNMED6,E_F1]= Nelder_Mead_ED6(E_HNMED6(i),S_HNMED6(i));
        
        %HNMED5
        %HNMED6
        %Vectors(j,1,1:Ns(i)+2)=HNMED5;
        %Vectors(j,1,1:Ns(i)+2)=HNMED6;
        %Vectors(j,1,1:Ns(i)+2)=C_LSHADE;
        %Solutions(j,:)=[HNMED5(Ns(i)+1)]
        %Solutions(j,:)=[HNMED6]
        Solutions(j,:)=[HJ(Ns(i)+1)]

        else
        %___________________Micro-Red______________________________________
        % El problema de la microred electrica se implemta diferente llamar
        % aqui los métodos para el problema
        [HNMED6,E_F1]= MR_HNMED6 (E_HNMED6(i),S_HNMED6(i),j);
        %[C_LSHADE,E_F1]= MR_CLSHADE(j);
         %C_LSHADE
        Solutions(j,:)=[HNMED6]
        %Solutions(j,:)=[C_LSHADE]
        end
        % Arreglos con tasas de convergencia para cada ejecución
%        E_F1s{j}=E_F1(1:end,:,:);
    end
    
    mejor=min(Solutions)
    peor=max(Solutions)
    mediana=median(Solutions)
    promedio=mean(Solutions)
    desvacion_estandar=std(Solutions)

    
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
    
    %*************** SALVANDO  ****************************
    % Grafica funciones de convergencian Experimento.
   % Name = ['/F-Datos_graficas_convergencia_CLSHADE_P_' num2str(i) '.xls'];
    Name = ['/E-Datos_graficas_convergencia_HNMED5_P_' num2str(i) '.xls'];
    %Name = ['/F-Datos_graficas_convergencia_HNMED6_P_' num2str(i) '.xls'];
    
    Nameexcel = [PathFolder Name];
    E_F1=squeeze(E_F1);
    xlswrite(Nameexcel,E_F1,1,'A1');
  
%     % Salvando resultados de Ejecuciones
%     Name = ['/F-Resultados_CLSHADE_Estadistica_Problema_' num2str(i) '.xls'];
%     %Name = ['/F-Resultados_HNMED5_Estadistica_Problema_' num2str(i) '.xls'];
%     %Name = ['/F-Resultados_HNMED6_Estadistica_Problema_' num2str(i) '.xls'];
%     Nameexcel = [PathFolder Name];
%     encabezado = {'CLSHADE'};
%     %encabezado = {'HNMED5'};
%     %encabezado = {'HNMED6'};
%     xlswrite(Nameexcel, encabezado, 1, 'B1');
%     xlswrite(Nameexcel, Solutions, 1, 'B2');
%     xlswrite(Nameexcel, {'Mejor'}, 1, 'A33');
%     xlswrite(Nameexcel, mejor, 1, 'B33');
%     xlswrite(Nameexcel, {'Peor'}, 1, 'A34');
%     xlswrite(Nameexcel, peor, 1, 'B34');
%     xlswrite(Nameexcel, {'Mediana'}, 1, 'A35');
%     xlswrite(Nameexcel, mediana, 1, 'B35');
%     xlswrite(Nameexcel, {'Promedio'}, 1, 'A36');
%     xlswrite(Nameexcel,promedio, 1, 'B36');
%     xlswrite(Nameexcel, {'Desv. Stan'}, 1, 'A37');
%     xlswrite(Nameexcel, desvacion_estandar, 1, 'B37');
%     xlswrite(Nameexcel, {'Mejor indivudio'}, 1, 'A38'); 
%     if i~=6
%     best=squeeze(Vectors(min_index(1),1,1:Ns(i)+2));
%     xlswrite(Nameexcel,best, 1, ['A39:A' num2str(40+Ns(i)+1)] ); 
%     else
%    % Eb el caso de la micro-Red solo se guarda el indice ya que los
%    % resultados de las ejecuciones se guardan a parte.
%     xlswrite(Nameexcel,min_index(1), 1, 'A39' );  
%     end
%     beep;
end

end