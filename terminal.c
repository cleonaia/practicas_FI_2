#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 512

float Mat[N][N];
float MatDD[N][N];
float V1[N], V2[N], V3[N], V4[N];

/* Inicialització */
void InitData(){
    int i,j;
    srand(221122);
    for( i = 0; i < N; i++ )
        for( j = 0; j < N; j++ ){
            Mat[i][j]=(((i*j)%3)?-1:1)(100.0f(rand()/(1.0f*RAND_MAX)));
            if ( (abs(i - j) <= 3) && (i != j))
                MatDD[i][j] = (((i*j)%3) ? -1 : 1)*(rand()/(1.0f*RAND_MAX));
            else if ( i == j )
                MatDD[i][j]=(((i*j)%3)?-1:1)(10000.0f(rand()/(1.0f*RAND_MAX)));
            else MatDD[i][j] = 0.0f;
        }
    for( i = 0; i < N; i++ ){
        V1[i]=(i<N/2)?(((i*j)%3)?-1:1)(100.0f(rand()/(1.0f*RAND_MAX))):0.0f;
        V2[i]=(i>=N/2)?(((i*j)%3)?-1:1)(100.0f(rand()/(1.0f*RAND_MAX))):0.0f;
        V3[i]=(((i*j)%5)?-1:1)(100.0f(rand()/(1.0f*RAND_MAX)));
    }
}

/* Impressió neta */
void PrintVect( float vect[N], int from, int numel ){
    for(int i = from; i < from + numel && i < N; i++)
        printf("%f ", vect[i]);
    printf("\n");
}
void PrintRow( float mat[N][N], int row, int from, int numel ){
    for(int i = from; i < from + numel && i < N; i++)
        printf("%f ", mat[row][i]);
    printf("\n");
}

/* Operacions vectorials/matricials */
void MultEscalar( float vect[N], float vectres[N], float alfa ){
    for(int i=0;i<N;i++) vectres[i] = vect[i]*alfa;
}
float Scalar( float vect1[N], float vect2[N] ){
    float s = 0.0f;
    for(int i=0;i<N;i++) s += vect1[i]*vect2[i];
    return s;
}
float Magnitude( float vect[N] ){
    return sqrtf(Scalar(vect,vect));
}
int Ortogonal( float vect1[N], float vect2[N] ){
    return (fabsf(Scalar(vect1,vect2)) < 1e-6f) ? 1 : 0;
}
/* Projecció segons enunciat: (u·v / |v|) * v */
void Projection( float vect1[N], float vect2[N], float vectres[N] ){
    float dot = Scalar(vect1,vect2);
    float mag = Magnitude(vect2);
    if (mag < 1e-12f){
        for(int i=0;i<N;i++) vectres[i]=0.0f;
        return;
    }
    float factor = dot / mag;
    for(int i=0;i<N;i++)
        vectres[i] = factor * vect2[i];
}

/* Normes de matriu */
float Infininorm( float M[N][N] ){
    float maxv = 0.0f;
    for(int i=0;i<N;i++){
        float sum=0.0f;
        for(int j=0;j<N;j++) sum += fabsf(M[i][j]);
        if(sum > maxv) maxv = sum;
    }
    return maxv;
}
float NormOne( float M[N][N] ){
    float maxc = 0.0f;
    for(int j=0;j<N;j++){
        float sum=0.0f;
        for(int i=0;i<N;i++) sum += fabsf(M[i][j]);
        if(sum > maxc) maxc = sum;
    }
    return maxc;
}
/* Frobenius en float */
float NormFrobenius( float M[N][N] ){
    float acc = 0.0f;
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            acc += M[i][j]*M[i][j];
    return sqrtf(acc);
}

/* Diagonal dominant */
int DiagonalDom( float M[N][N] ){
    for(int i=0;i<N;i++){
        float diag = fabsf(M[i][i]);
        float sum = 0.0f;
        for(int j=0;j<N;j++)
            if(j!=i) sum += fabsf(M[i][j]);
        if(diag < sum) return 0;
    }
    return 1;
}

/* Matriu * vector */
void Matriu_x_Vector( float M[N][N], float vect[N], float vectres[N] ){
    for(int i=0;i<N;i++){
        float acc=0.0f;
        for(int j=0;j<N;j++)
            acc += M[i][j]*vect[j];
        vectres[i]=acc;
    }
}

/* Jacobi bàsic (una iteració a partir de x_in) */
int JacobiOneIteration(float M[N][N], float b[N], float x_in[N], float x_out[N]){
    if(!DiagonalDom(M)) return 0;
    for(int i=0;i<N;i++){
        float diag = M[i][i];
        if (fabsf(diag) < 1e-12f) return 0;
        float sum=0.0f;
        for(int j=0;j<N;j++)
            if(j!=i) sum += M[i][j]*x_in[j];
        x_out[i] = (b[i]-sum)/diag;
    }
    return 1;
}

/* Jacobi (iters) des de x=0 */
int JacobiIter(float M[N][N], float b[N], unsigned iters, float x_out[N]){
    if(!DiagonalDom(M)) return 0;
    float x_curr[N]={0.0f};
    float x_next[N];
    for(unsigned k=0;k<iters;k++){
        if(!JacobiOneIteration(M,b,x_curr,x_next)) return 0;
        for(int i=0;i<N;i++) x_curr[i]=x_next[i];
    }
    for(int i=0;i<N;i++) x_out[i]=x_curr[i];
    return 1;
}

/* Error com ||x(k+1) − x(k)||_2 (L2 global a tots els N elements) */
float DiffL2_all(const float a[N], const float b[N]){
    float acc = 0.0f;
    for (int i = 0; i < N; i++){
        float d = a[i] - b[i];
        acc += d * d;
    }
    return sqrtf(acc);
}

/* Punt Extra */
void ResidualVec(float M[N][N], float x[N], float b[N], float r[N]){
    for(int i=0;i<N;i++){
        float acc=0.0f;
        for(int j=0;j<N;j++) acc += M[i][j]*x[j];
        r[i] = acc - b[i];
    }
}
float Vector2Norm(float v[N]){
    double acc=0.0;
    for(int i=0;i<N;i++) acc += (double)v[i]*(double)v[i];
    return (float)sqrt(acc);
}
float JacobiThetaInf(float A[N][N]){
    float theta = 0.0f;
    for(int i=0;i<N;i++){
        float diag = fabsf(A[i][i]);
        if (diag < 1e-20f) return 1.0f;
        float rowSum = 0.0f;
        for(int j=0;j<N;j++) if (j!=i) rowSum += fabsf(A[i][j]);
        float val = rowSum / diag;
        if (val > theta) theta = val;
    }
    if (theta >= 1.0f) theta = 0.999999f; /* evitar div. per zero a la cota */
    return theta;
}

int main(void){
    InitData();

    /* A */
    printf("V1 del 0 al 9 i del 256 al 265:\n");
    PrintVect(V1,0,10);
    PrintVect(V1,256,10);
    printf("\n");

    printf("V2 del 0 al 9 i del 256 al 265:\n");
    PrintVect(V2,0,10);
    PrintVect(V2,256,10);
    printf("\n");

    printf("V3 del 0 al 9 i del 256 al 265:\n");
    PrintVect(V3,0,10);
    PrintVect(V3,256,10);
    printf("\n");

    /* B */
    printf("Mat fila 0 i fila 100 del 0 al 9:\n");
    PrintRow(Mat,0,0,10);
    PrintRow(Mat,100,0,10);
    printf("\n");

    /* C */
    printf("MatDD fila 0 del 0 al 9 i fila 100 del 95 al 104:\n");
    PrintRow(MatDD,0,0,10);
    PrintRow(MatDD,100,95,10);
    printf("\n");

    /* Normes i diagonalitat */
    printf("Infininorma de Mat = %.3f\n", Infininorm(Mat));
    printf("Norma ú de Mat = %.3f\n", NormOne(Mat));
    printf("Norma de Frobenius de Mat = %.3f\n", NormFrobenius(Mat));
    printf("La matriu Mat %s diagonal dominant\n", DiagonalDom(Mat)? "és" : "no és");
    printf("\n");

    printf("Infininorma de MatDD = %.3f\n", Infininorm(MatDD));
    printf("Norma ú de MatDD = %.3f\n", NormOne(MatDD));
    printf("Norma de Frobenius de MatDD = %.3f\n", NormFrobenius(MatDD));
    printf("La matriu MatDD %s diagonal dominant\n", DiagonalDom(MatDD)? "és" : "no és");
    printf("\n");

    /* Escalars, magnituds i ortogonalitat */
    printf("Escalar <V1,V2> = %f Escalar <V1,V3> = %f Escalar <V2,V3> = %f\n",
           Scalar(V1,V2), Scalar(V1,V3), Scalar(V2,V3));
    printf("Magnitud V1,V2 i V3 = %f %f %f\n",
           Magnitude(V1), Magnitude(V2), Magnitude(V3));
    printf("V1 i V2 %s ortogonals\n", Ortogonal(V1,V2)? "són" : "no són");
    printf("\n");

    /* V3 * 2 */
    float V3x2[N];
    MultEscalar(V3,V3x2,2.0f);
    printf("Els elements 0 al 9 i 256 al 265 del resultat de multiplicar V3x2.0 són:\n");
    PrintVect(V3x2,0,10);
    PrintVect(V3x2,256,10);
    printf("\n");

    /* Projeccions */
    float ProjV2onV3[N], ProjV1onV2[N];
    Projection(V2,V3,ProjV2onV3);
    Projection(V1,V2,ProjV1onV2);
    printf("Els elements 0 a 9 del resultat de la projecció de V2 sobre V3 són:\n");
    PrintVect(ProjV2onV3,0,10);
    printf("Els elements 0 a 9 del resultat de la projecció de V1 sobre V2 són:\n");
    PrintVect(ProjV1onV2,0,10);
    printf("\n");

    /* Mat * V2 */
    float MatV2[N];
    Matriu_x_Vector(Mat,V2,MatV2);
    printf("Els elements 0 a 9 del resultat de la multiplicació de Mat per v2 són:\n");
    PrintVect(MatV2,0,10);
    printf("\n");

    /* Jacobi sobre MatDD*/
    if (!DiagonalDom(MatDD)){
        printf("La matriu MatDD no és diagonal dominant, no es pot aplicar Jacobi\n");
    } else {
        float x1[N], x2[N], x1000[N], x1001[N];

        /* x^(1) i x^(2) */
        JacobiIter(MatDD, V3, 1, x1);
        JacobiOneIteration(MatDD, V3, x1, x2);

        /* x^(1000) i x^(1001) */
        JacobiIter(MatDD, V3, 1000, x1000);
        JacobiOneIteration(MatDD, V3, x1000, x1001);

    /* Error com ||x(k+1) − x(k)||_2 (L2 global) */
    float error1    = DiffL2_all(x2,     x1);
    float error1000 = DiffL2_all(x1001,  x1000);

        printf("Els elements 0 a 9 de la solució (1 iter) del sistema d'equacions són:\n");
        PrintVect(x1,0,10);
    printf("Error %.6f\n", error1);
        printf("\n");

        printf("Els elements 0 a 9 de la solució (1000 iters) del sistema d'equacions són:\n");
        PrintVect(x1000,0,10);
    printf("Error %.6f\n", error1000);
        printf("\n");

        /* —————————————————— Punt extra —————————————————— */
        float r[N];
        float b2 = Magnitude(V3);

        /* Residus i relatius */
        ResidualVec(MatDD, x1,    V3, r);
        float r2_1     = Vector2Norm(r);
        float rel2_1   = (b2>0.0f)? (r2_1/b2) : 0.0f;

        ResidualVec(MatDD, x1000, V3, r);
        float r2_1000   = Vector2Norm(r);
        float rel2_1000 = (b2>0.0f)? (r2_1000/b2) : 0.0f;

        /* Cota de Jacobi amb ||G||_inf */
        float theta = JacobiThetaInf(MatDD);
        float bound1     = (theta/(1.0f - theta)) * error1;
        float bound1000  = (theta/(1.0f - theta)) * error1000;

        /* Punt extra */
        printf("\n======================== Punt extra — qualitat Jacobi ========================\n");
        printf("||G||_inf = %.6f\n", theta);
        printf("  (1 iter)    ||r||2 = %.6f    rel(||r||2) = %.6f    cota_inf(error) ≤ %.6f\n",
               r2_1,    rel2_1,    bound1);
        printf("  (1000 it.)  ||r||2 = %.6f    rel(||r||2) = %.6f    cota_inf(error) ≤ %.6f\n",
               r2_1000, rel2_1000, bound1000);
        printf("=============================================================================\n\n");
    }

    /* Jacobi sobre Mat */
    if (!DiagonalDom(Mat)){
        printf("La matriu M no és diagonal dominant, no es pot aplicar Jacobi\n");
    } else {
        float dummy[N];
        JacobiIter(Mat, V3, 10, dummy);
        PrintVect(dummy,0,10);
    }

    return 0;
}

/*
Punt extra:

•⁠  ⁠Mesurem la qualitat de la solució x d’Ax=b de dues maneres:
  1) Residual: r = A x − b. Mostrem ||r||2 i el residual relatiu ||r||2 / ||b||2.
  2) Cota a posteriori de Jacobi amb ∞-norma:
     G = −D^{-1}(L+U), θ = ||G||inf = max_i Σ{j≠i} |a_ij|/|a_ii|.
     Si θ < 1, llavors ||x* − x(k)||_inf ≤ (θ/(1−θ)) · ||x(k+1) − x(k)||_inf.
•⁠  ⁠L’error que es mostra per a l’apartat J s’ha calculat com la norma L2 nomésdels 10 primers elements.
•⁠  ⁠Hem fet una complicació diferent per evitar diferències:
  gcc main.c -O0 -lm
*/
