#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 512
#define EPS_ORTO 1e-6f
#define EPS_DIAG 1e-8f

static float Mat[N][N], MatDD[N][N];
static float V1[N], V2[N], V3[N];

static inline float randf(float scale){ return scale * (rand() / (float)RAND_MAX); }

void InitData(void){
    srand(221122);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            Mat[i][j] = (((i*j)%3)?-1.0f:1.0f)*randf(100.0f);
            if(i==j){
                MatDD[i][j] = (((i*j)%3)?-1.0f:1.0f)*randf(10000.0f);
            } else if (abs(i-j)<=3){
                MatDD[i][j] = (((i*j)%3)?-1.0f:1.0f)*randf(1.0f);
            } else {
                MatDD[i][j] = 0.0f;
            }
        }
        V1[i] = (i < N/2) ? (((i)%3)?-1.0f:1.0f)*randf(100.0f) : 0.0f;
        V2[i] = (i >= N/2) ? (((i)%3)?-1.0f:1.0f)*randf(100.0f) : 0.0f;
        V3[i] = (((i)%5)?-1.0f:1.0f)*randf(100.0f);
    }
}

void PrintVect(const float v[N], int from, int num){
    for(int i=from;i<from+num && i<N;i++)
        printf("Vect[%d]= %f\n", i, v[i]);
}

void PrintRow(const float M[N][N], int row, int from, int num){
    printf("Row %d: ", row);
    for(int j=from;j<from+num && j<N;j++)
        printf("%f ", M[row][j]);
    putchar('\n');
}

void MultEscalar(const float v[N], float vr[N], float a){
    for(int i=0;i<N;i++) vr[i]=v[i]*a;
}

float Scalar(const float a[N], const float b[N]){
    float s=0.0f;
    for(int i=0;i<N;i++) s += a[i]*b[i];
    return s;
}

float Magnitude(const float v[N]){
    return sqrtf(Scalar(v,v));
}

int Ortogonal(const float a[N], const float b[N]){
    return fabsf(Scalar(a,b)) < EPS_ORTO;
}

void Projection(const float u[N], const float v[N], float pr[N]){
    float denom = Scalar(v,v);
    float coef  = denom>0.0f ? Scalar(u,v)/denom : 0.0f;
    for(int i=0;i<N;i++) pr[i]=coef*v[i];
}

float Infininorm(const float M[N][N]){
    float mx=0.0f;
    for(int i=0;i<N;i++){
        float s=0.0f;
        for(int j=0;j<N;j++) s += fabsf(M[i][j]);
        if(s>mx) mx=s;
    }
    return mx;
}

float NormFrobenius(const float M[N][N]){
    double acc=0.0; // usar doble para menor error acumulado
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            acc += (double)M[i][j] * (double)M[i][j];
    return (float)sqrt(acc);
}

int DiagonalDom(const float M[N][N]){
    for(int i=0;i<N;i++){
        float diag=fabsf(M[i][i]), off=0.0f;
        for(int j=0;j<N;j++) if(j!=i) off += fabsf(M[i][j]);
        if(diag < off) return 0;
    }
    return 1;
}

void Matriu_x_Vector(const float M[N][N], const float v[N], float r[N]){
    for(int i=0;i<N;i++){
        float s=0.0f;
        for(int j=0;j<N;j++) s += M[i][j]*v[j];
        r[i]=s;
    }
}

float ResidualNorm(const float M[N][N], const float x[N], const float b[N]){
    double acc=0.0;
    for(int i=0;i<N;i++){
        double row=0.0;
        for(int j=0;j<N;j++) row += M[i][j]*x[j];
        double diff = row - b[i];
        acc += diff*diff;
    }
    return (float)sqrt(acc);
}

int Jacobi(const float M[N][N], const float b[N], float x_out[N], unsigned iters){
    if(!DiagonalDom(M)) return 0;
    float x[N]={0};
    for(unsigned k=0;k<iters;k++){
        float next[N];
        for(int i=0;i<N;i++){
            if(fabsf(M[i][i])<EPS_DIAG) return 0;
            float s=0.0f;
            for(int j=0;j<N;j++) if(j!=i) s+= M[i][j]*x[j];
            next[i] = (b[i]-s)/M[i][i];
        }
        for(int i=0;i<N;i++) x[i]=next[i];
    }
    for(int i=0;i<N;i++) x_out[i]=x[i];
    return 1;
}

int main(void){
    InitData();

    puts("V1 (0-9) y (256-265):");
    PrintVect(V1,0,10); PrintVect(V1,256,10);
    puts("V2 (0-9) y (256-265):");
    PrintVect(V2,0,10); PrintVect(V2,256,10);
    puts("V3 (0-9) y (256-265):");
    PrintVect(V3,0,10); PrintVect(V3,256,10);

    puts("Mat filas 0 y 100 (0-9):");
    PrintRow(Mat,0,0,10); PrintRow(Mat,100,0,10);

    puts("MatDD fila 0 (0-9) y fila 100 (95-104):");
    PrintRow(MatDD,0,0,10); PrintRow(MatDD,100,95,10);

    printf("Inf-norm Mat = %.3f\n", Infininorm(Mat));
    printf("Frobenius Mat = %.3f\n", NormFrobenius(Mat));
    printf("Mat %s diagonal dominante\n", DiagonalDom(Mat)?"es":"no es");

    printf("Inf-norm MatDD = %.3f\n", Infininorm(MatDD));
    printf("Frobenius MatDD = %.3f\n", NormFrobenius(MatDD));
    printf("MatDD %s diagonal dominante\n", DiagonalDom(MatDD)?"es":"no es");

    printf("<V1,V2>= %f\n", Scalar(V1,V2));
    printf("<V1,V3>= %f\n", Scalar(V1,V3));
    printf("<V2,V3>= %f\n", Scalar(V2,V3));

    printf("||V1|| ||V2|| ||V3|| = %f %f %f\n", Magnitude(V1), Magnitude(V2), Magnitude(V3));

    printf("V1 y V2 %sson ortogonales\n", Ortogonal(V1,V2)?"":"no ");
    printf("V1 y V3 %sson ortogonales\n", Ortogonal(V1,V3)?"":"no ");
    printf("V2 y V3 %sson ortogonales\n", Ortogonal(V2,V3)?"":"no ");

    float V3x2[N]; MultEscalar(V3,V3x2,2.0f);
    puts("V3*2 (0-9) y (256-265):");
    PrintVect(V3x2,0,10); PrintVect(V3x2,256,10);

    float ProjV2onV3[N], ProjV1onV2[N];
    Projection(V2,V3,ProjV2onV3);
    Projection(V1,V2,ProjV1onV2);
    puts("Proy V2 sobre V3 (0-9):"); PrintVect(ProjV2onV3,0,10);
    puts("Proy V1 sobre V2 (0-9):"); PrintVect(ProjV1onV2,0,10);

    float MatV2[N]; Matriu_x_Vector(Mat,V2,MatV2);
    puts("Mat*V2 (0-9):"); PrintVect(MatV2,0,10);

    float sol1[N], sol1000[N];
    if(Jacobi(MatDD,V3,sol1,1)){
        puts("Jacobi MatDD (1 iter) (0-9):"); PrintVect(sol1,0,10);
        printf("Residual 1 iter = %f\n", ResidualNorm(MatDD,sol1,V3));
    }
    if(Jacobi(MatDD,V3,sol1000,1000)){
        puts("Jacobi MatDD (1000 iters) (0-9):"); PrintVect(sol1000,0,10);
        printf("Residual 1000 iters = %f\n", ResidualNorm(MatDD,sol1000,V3));
    }

    float solMat[N];
    if(Jacobi(Mat,V3,solMat,1000)){
        puts("Jacobi Mat (1000 iters) (0-9):"); PrintVect(solMat,0,10);
    } else {
        puts("Mat no es diagonal dominante; Jacobi no aplicable.");
    }

    return 0;
}

