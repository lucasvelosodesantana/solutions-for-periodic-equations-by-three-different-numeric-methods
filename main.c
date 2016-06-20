#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

   typedef struct var {
       double a0, a1, a2, a3, a4, a5;
       int Nmax;
       double e1, e2;
       double a, b;
       } Variaveis;

double fx(double x, double a0, double a1, double a2, double a3, double a4, double a5);
double fxd1(double x, double a0, double a1, double a2, double a3, double a4, double a5);
double fxd2(double x, double a0, double a1, double a2, double a3, double a4, double a5);
void newton(double a0, double a1, double a2, double a3, double a4, double a5, int Nmax, double e1, double e2, double a, double b, int n);
void bissecao(double a0, double a1, double a2, double a3, double a4, double a5, int Nmax, double e1, double e2, double a, double b);
void halley(double a0, double a1, double a2, double a3, double a4, double a5, int Nmax, double e1, double e2, double a, double b);

int main(){
    Variaveis t1;
    FILE* ArqEntrada;
    ArqEntrada = fopen("dados.txt", "r");
    int lidos, n = 1;
    if(ArqEntrada == NULL){
        printf("Erro ao abrir o arquivo!\n");
    }
    else{
        while(!feof(ArqEntrada)){
            lidos = fscanf(ArqEntrada, "%lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf", &t1.a0, &t1.a1, &t1.a2, &t1.a3, &t1.a4, &t1.a5, &t1.Nmax, &t1.e1, &t1.e2, &t1.a, &t1.b);
            if(lidos > 0){
                newton(t1.a0 ,t1.a1, t1.a2, t1.a3, t1.a4, t1.a5,t1.Nmax, t1.e1, t1.e2, t1.a, t1.b, n);
                halley(t1.a0 ,t1.a1, t1.a2, t1.a3, t1.a4, t1.a5,t1.Nmax, t1.e1, t1.e2, t1.a, t1.b);
                bissecao(t1.a0, t1.a1, t1.a2, t1.a3, t1.a4, t1.a5, t1.Nmax, t1.e1, t1.e2, t1.a, t1.b);
                n = n + 1;
            }
        }
    }
    fclose(ArqEntrada);
    return 0;
}
double fx(double x, double a0, double a1, double a2, double a3, double a4, double a5){
    double rec = a0*cos(a1*x) + a2*sin(a3*x) + exp(a4*x) + a5;
    return rec;
}
double fxd1(double x, double a0, double a1, double a2, double a3, double a4, double a5){
    double rec = -sin(a1*x)*a0*a1 + cos(a3*x)*a2*a3 + a4*exp(a4*x);
    return rec;
}
double fxd2(double x, double a0, double a1, double a2, double a3, double a4, double a5){
    double rec = -cos(a1*x)*a0*a1*a1 - sin(a3*x)*a2*a3*a3 + exp(a4*x)*a4*a4;
    return rec;
}

void newton(double a0, double a1, double a2, double a3, double a4, double a5, int Nmax, double e1, double e2, double a, double b, int n){
    FILE* ArqSaida;
    ArqSaida = fopen("resultados.txt", "a");
    int k = 0;
    double x0 = (a+b)/2, x1 = 0, fx0, fx0d1;
    fx0 = fx(x0, a0, a1, a2, a3, a4, a5);
    if(ArqSaida == NULL){
        printf("Erro ao criar o arquivo\n");
    }
    else{
        fprintf(ArqSaida, "Entrada: %d\nMétodo: Newton-Raphson\n\n", n);
    }
    if(fabs(fx0) > e2){
        while(((fabs(x0 - x1)>e1)||(fabs(fx0)>e2))&&(k < Nmax)){
            k = k + 1;
            x1 = x0;
            fx0d1 = fxd1(x0, a0, a1, a2, a3, a4, a5);
            x0 = x1 - (fx0)/(fx0d1);
            fx0 = fx(x0, a0, a1, a2, a3, a4, a5);
        }
        if(ArqSaida == NULL){
            printf("Erro ao criar o arquivo\n");
        }
        else{
            if(fabs(fx0) < e2){
                fprintf(ArqSaida, "Convergiu!\n");
            }
            else{
                fprintf(ArqSaida, "Não Convergiu!\n");
            }
            fprintf(ArqSaida, "Número de iterações: %d\nRaiz final: %.5lf\n", k, x0);
            fprintf(ArqSaida, "|x_(%d)-x_(%d)|: %.5lf\n", k, k-1, fabs(x0 - x1));
            fprintf(ArqSaida, "|f(x_%d)|: %.5lf\n\n", k, fabs(fx0));
        }
    }
    else{
        if(ArqSaida == NULL){
            printf("Erro ao criar o arquivo\n");
        }
        else{
            fprintf(ArqSaida, "Convergiu!\n");
            fprintf(ArqSaida, "Número de iterações: 0\nRaiz final: %.5lf\n", x0);
            fprintf(ArqSaida, "|x_(i)-x_(i)|: Obtido sem iterações\n");
            fprintf(ArqSaida, "|f(x_0)|: %.5lf\n\n", fx0);
        }
    }
    fclose(ArqSaida);
}

void bissecao(double a0, double a1, double a2, double a3, double a4, double a5, int Nmax, double e1, double e2, double a, double b){
    FILE* ArqSaida;
    ArqSaida = fopen("resultados.txt", "a");
    int k = 0;
    double x0 = (a+b)/2, x1 = 0, c = (b - a), fx0;
    fx0 = fx(x0, a0, a1, a2, a3, a4, a5);
    if(ArqSaida == NULL){
        printf("Erro ao criar o arquivo\n");
    }
    else{
        fprintf(ArqSaida, "Método: Bisseção\n\n");
    }
    if(fabs(fx0)> e2){
        while(((fabs(c)>e1)||(fabs(fx0)>e2))&&(k < Nmax)){
            k = k + 1;
            if((fx(a, a0, a1, a2, a3, a4, a5)*fx0) == 0){
                x1 = x0;
                b = x0;
            }
            if((fx(a, a0, a1, a2, a3, a4, a5)*fx0) < 0){
                x1 = x0;
                b = x0;
            }
            if(((fx(a, a0, a1, a2, a3, a4, a5))*fx0) > 0){
                x1 = x0;
                a = x0;
            }
            c = b - a;
            x0 = (a+b)/2;
            fx0 = fx(x0, a0, a1, a2, a3, a4, a5);
        }
        if(ArqSaida == NULL){
            printf("Erro ao criar o arquivo\n");
        }
        else{
            if(fabs(fx0) < e2){
                fprintf(ArqSaida, "Convergiu!\n");
            }
            else{
                fprintf(ArqSaida, "Não Convergiu!\n");
            }
            fprintf(ArqSaida, "Número de iterações: %d\nRaiz final: %.5lf\n", k, x0);
            fprintf(ArqSaida, "|x_(%d)-x_(%d)|: %.5lf\n", k, k-1, fabs(x0 - x1));
            fprintf(ArqSaida, "|f(x_%d)|: %.5lf\n\n", k, fabs(fx0));
        }
    }
    else{
        if(ArqSaida == NULL){
            printf("Erro ao criar o arquivo\n");
        }
        else{
            fprintf(ArqSaida, "Convergiu!\n");
            fprintf(ArqSaida, "Número de iterações: 0\nRaiz final: %.5lf\n", x0);
            fprintf(ArqSaida, "|x_(i)-x_(i)|: Obtido sem iterações\n");
            fprintf(ArqSaida, "|f(x_0)|: %.5lf\n\n", fx0);
        }
    }
    fclose(ArqSaida);
}

void halley(double a0, double a1, double a2, double a3, double a4, double a5, int Nmax, double e1, double e2, double a, double b){
    FILE* ArqSaida;
    ArqSaida = fopen("resultados.txt", "a");
    int k = 0;
    double x0 = (a+b)/2, x1 = 0, fx0, fx0d1, fx0d2;
    fx0 = fx(x0, a0, a1, a2, a3, a4, a5);
    if(ArqSaida == NULL){
        printf("Erro ao criar o arquivo\n");
    }
    else{
        fprintf(ArqSaida, "Método: Halley\n\n");
    }
    if(fabs(fx0) > e2){
        while(((fabs(x0 - x1)>e1)||(fabs(fx0)>e2))&&(k < Nmax)){
            k = k + 1;
            x1 = x0;
            fx0d1 = fxd1(x0, a0, a1, a2, a3, a4, a5);
            fx0d2 = fxd2(x0, a0, a1, a2, a3, a4, a5);
            x0 = x1 - (2*fx0*fx0d1)/((2*fx0d1*fx0d1)-(fx0*fx0d2));
            fx0 = fx(x0, a0, a1, a2, a3, a4, a5);
        }
        if(ArqSaida == NULL){
            printf("Erro ao criar o arquivo\n");
        }
        else{
            if(fabs(fx0) < e2){
                fprintf(ArqSaida, "Convergiu!\n");
            }
            else{
                fprintf(ArqSaida, "Não Convergiu!\n");
            }
            fprintf(ArqSaida, "Número de iterações: %d\nRaiz final: %.5lf\n", k, x0);
            fprintf(ArqSaida, "|x_(%d)-x_(%d)|: %.5lf\n", k, k-1, fabs(x0 - x1));
            fprintf(ArqSaida, "|f(x_%d)|: %.5lf\n\n", k, fabs(fx0));
        }
    }
    else{
        if(ArqSaida == NULL){
            printf("Erro ao criar o arquivo\n");
        }
        else{
            fprintf(ArqSaida, "Convergiu!\n");
            fprintf(ArqSaida, "Número de iterações: 0\nRaiz final: %.5lf\n", x0);
            fprintf(ArqSaida, "|x_(i)-x_(i)|: Obtido sem iterações\n");
            fprintf(ArqSaida, "|f(x_0)|: %.5lf\n\n", fx0);
        }
    }
    fclose(ArqSaida);
}
