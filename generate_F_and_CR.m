% genera el factor de escala F y la probabilidad de mutación CR
function[F,CR]=generate_F_and_CR()
    F=rand();
    while F<0.3 || F>0.9
        F=rand();
    end
                  
    CR=rand();
    while CR<0.8
        CR=rand();
    end
end