function[cost,cvs]= MR_cost(P)
    P1=P(1);
    P2=P(2);
    P3=P(3);
    P4=P(4);
    global SOC_0;
    [sum_P1,sum_P2,sum_P3,sum_P4]=get_sum_P1_P2_P2_P4();
    t=get_hour();
    [w1, w2, w3, w4] = deal(0.25);
    cf=1;
    cost=0;
    cvs=0;
    %%parámetros bateria para restricción SOC
    %n_B=0.85; % Round trip efficiency
    n_c=0.85;
    n_d=1;
    Bc_max=2000; %Capacidad máxima
    alpha_c=n_c/Bc_max;
    alpha_d=n_d/Bc_max;
    SOC_max=0.95;    %Estado de carga mínimo de la bateria
    SOC_min=0.4;    %Estado de carga máximo de la bateria
   
    %*****Función objetivo******* 
    cost=cost+cf*w1*F_1(P1)+w2*F_2(P2)-w3*F_3(P3)-w4*F_4(P4); 
    %Balance de potencia: 
    b_NNO=-SOC_min+SOC_0-alpha_d*sum_P2+alpha_c*sum_P3+alpha_c*sum_P4;
    b_DO=SOC_max-SOC_0+alpha_d*sum_P2-alpha_c*sum_P3-alpha_c*sum_P4; 
    %La suma de las potencias suministradas por los generadores 
    % debe ser igual a la carga total del sistema, esto es:
     constrain=abs(P1 + P2 + P3 + P4-get_PL());
     %if(constrain>0.000001)
     cvs=cvs+constrain;
    % end
     Lizq1=alpha_d*P2-alpha_c*P3-alpha_c*P4;
     cvs=cvs+Lizq1-b_NNO;
     Lizq2=-alpha_d*P2+alpha_c*P3+alpha_c*P4;
     cvs=cvs+Lizq2-b_DO;
    %end
end

function [C_P1]=F_1(P1)
    alfa= 14.88; 
    beta=0.3;
    gamma= 0.000435;
    C_P1 = alfa + beta*P1 + gamma*(P1^2);
end
function [C_P2]=F_2(P2)
    C_P2=119*P2;
end
function [C_P3]=F_3(P3)
    C_P3=545.016*P3 ;
end
function [C_P4]=F_4(P4)
   C_P4= 152.616*P4; 
end
