#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void methodKramera(double* x_alpha, double alpha, double** Ah, double* f);
double residual(double* x_alpha, double alpha, double resid);
void Init(double** Ah, double* f, double* res);

int main(){
    double* x;
    double** Ah;
    double* f;
    double a = 0.0, b = 1.0, alpha, resid;
    double acX_1 = 0.48, acX_2 = 0.36;
    double infelicity;
    f = (double*)malloc(sizeof(double)*2);
    x = (double*)malloc(sizeof(double)*2);
    Ah = (double**)malloc((sizeof(double*)*2));
    for(int i = 0; i < 2; ++i){
        Ah[i] = (double*)malloc(sizeof (double)*2);
    }

    if(x == NULL || f == NULL || Ah == NULL){
        printf("error memory\n");
        return -1;
    }

    Init(Ah, f, &resid);

    while (fabs(b-a) > 0.001) {
        double resB, resAl;
        methodKramera(x, b, Ah, f);
        resB = residual(x, b, resid);
        alpha = (a + b)/2;
        methodKramera(x, alpha, Ah, f);
        resAl = residual(x, alpha, resid);
        if(resB*resAl < 0.0){
            a = alpha;
        }else if(resB*resAl > 0.0){
            b = alpha;
        }else{
            if(resB == 0.0){
                alpha = b;
                break;
            }else{
                break;
            }
        }
    }

    alpha = (a+b)/2;
    methodKramera(x, alpha, Ah, f);
    residual(x, alpha, resid);
    printf("\n| Точное решение, найденное вручную |\n");
    printf("|-----------------------------------|\n");
    printf("| x[0] = %lf | x[1] = %lf |\n\n", acX_1, acX_2);
    printf("| Приближённое решение, найденное с помощью метода Лаврентьева |\n");
    printf("|--------------------------------------------------------------|\n");
    printf("| x[0] = %lf | x[1] = %lf | альфа: %lf\t       |\n\n", x[0], x[1], alpha);
    infelicity = fabs(acX_1 - x[0]);
    printf("Погрешность первого приближённого корня: %lf\n", infelicity);
    infelicity = fabs(acX_2 - x[1]);
    printf("Погрешность второго приближённого корня: %lf\n\n", infelicity);
    free(x);
    free(f);
    for(int i = 0; i < 2; ++i){
        free(Ah[i]);
    }
    free(Ah);

    return 0;
}

void methodKramera(double* x_alpha, double alpha, double** Ah, double* f){
    double delta_1, delta_2, delta_3;

    Ah[0][0] = 4.0 + 0.001  + alpha;
    Ah[1][1] = 2.25 + 0.003 + alpha;

    delta_1 = (Ah[0][0])*(Ah[1][1]) - Ah[0][1]*(Ah[1][0]);
    if(delta_1 == 0.0){
        printf("error, division on zero\n");
        return;
    }
    delta_2 = f[0]*(Ah[1][1]) - Ah[0][1]*(f[1]);
    delta_3 = (Ah[0][0])*(f[1]) - f[0]*(Ah[1][0]);

    x_alpha[0] = delta_2/delta_1;
    x_alpha[1] = delta_3/delta_1;
}

double residual(double* x_alpha, double alpha, double resid){

    double norma = sqrt(pow(x_alpha[0], 2) + pow(x_alpha[1], 2));
    double res = alpha*norma - resid;
    return res;
}

void Init(double** Ah, double* f, double* res){
    double K = 3.0, c = 0.7;
    double h = 0.005, delta = 0.003;
    double p = (1.0/2.0);
    double h_1 = 0.001, h_2 = 0.003, h_3 = -0.003;
    double deltaF_1 = -0.002, deltaF_2 = 0.002;

    *res = K*pow(delta+c*h, p);
    f[0] = 3.0 + deltaF_1;
    f[1] = 2.25 + deltaF_2;

    Ah[0][0] = 4.0 + h_1;
    Ah[1][1] = 2.25 + h_2;
    Ah[0][1] = 3.0 + h_3;
    Ah[1][0] = 3.0 + h_3;
}
