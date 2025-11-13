#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 512

float Mat[N][N];
float MatDD[N][N];
float V1[N], V2[N], V3[N], V4[N];

void InitData(){
int i,j;
srand(221122);
for( i = 0; i < N; i++ )
for( j = 0; j < N; j++ ){
Mat[i][j]=(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
if ( (abs(i - j) <= 3) && (i != j))
MatDD[i][j] = (((i*j)%3) ? -1 : 1)*(rand()/(1.0*RAND_MAX));
else if ( i == j )
MatDD[i][j]=(((i*j)%3)?-1:1)*(10000.0*(rand()/(1.0*RAND_MAX)));
else MatDD[i][j] = 0.0;
}
for( i = 0; i < N; i++ ){
V1[i]=(i<N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
V2[i]=(i>=N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
V3[i]=(((i*j)%5)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
}
}

void PrintVect( float vect[N], int from, int numel )
{
// 1. Mostrar un cert nombre d'elements d'un vector a partir d'una posicio donada
    for (int i = from; i < from + numel && i < N; i++ )
        printf( "Vect[%d]= %f\n", i, vect[i] );
}

// 2. Mostrar un cert nombre d'elements d'una fila d'una matriu a partir d'una posicio donada
void PrintRow( float mat[N][N], int row, int from, int numel )
{
    
    printf("Row %d: ", row);
    for(int i = from; i < from + numel && i < N; i++ )
        printf("%f ", mat[row][i]);
    printf("\n");
}

// 3. Multiplicar un vector per un escalar
void MultEscalar( float vect[N], float vectres[N], float alfa )
{
    for(int i = 0; i < N; i++ )
        vectres[i] = vect[i] * alfa;
}

// 4. Calcular el producte escalar de dos vectors
float Scalar( float vect1[N], float vect2[N] )
{
    float dotProduct = 0;
    for (int i = 0; i < N; i++) {
        dotProduct += vect1[i] * vect2[i];
    }
    return dotProduct;

}

float Magnitude(float vect[N] ){
    float sum = 0.0;
    for(int i = 0; i < N; i++){
        sum += vect[i] * vect[i];
    }
    return sqrt(sum);
}

int Ortogonal(float vect1[N], float vect2[N]){
    float dot_product = 0.0;
    for(int i = 0; i < N; i++){
        dot_product += vect1[i] * vect2[i];
    }
    return (fabs(dot_product) < 1e-6) ? 1 : 0;
}

void Projection(float vect1[N], float vect2[N], float vectres[N]){
    float dot_product = 0.0;
    float norm2 = 0.0; 
    
    for(int i = 0; i < N; i++){
        dot_product += vect1[i] * vect2[i];
        norm2 += vect2[i] * vect2[i];
    }
    float project_scalar = dot_product / norm2;
    for(int i = 0; i < N; i++){
        vectres[i] = project_scalar * vect2[i];
    }
}

float Infininorm(float M[N][N]){
    float max_val = 0.0;

    for(int i = 0; i < N; i++){
        float sum_row = 0.0;
        for(int j = 0; j < N; j++){
            sum_row += fabs(M[i][j]);
        }
        if(sum_row > max_val){
            max_val = sum_row;
        }
    }
    return max_val;
}

float NormFrobenius(float M[N][N]){
    float sum_squares = 0.0;
    for(int i = 0; i < N; i ++){
        for(int j = 0; j < N; j++){
            sum_squares += M[i][j] * M[i][j];
        }
    }
    return sqrt(sum_squares);
}

int DiagonalDom(float M[N][N]){
    for(int i = 0; i < N; i++){
        float diag_elem = fabs(M[i][i]);
        float sum_row = 0.0;
        for(int j = 0; j < N; j++){
            if(j != i){
                sum_row += fabs(M[i][j]);
            }
        }
        if(diag_elem < sum_row){
            return 0; // Diagonal no dominant
        }
    }
    return 1; // Diagonal dominant
}

void Matriu_x_Vector( float M[N][N], float vect[N], float vectres[N] ) {
    for(int i = 0; i < N; i++){
        vectres[i] = 0.0f;
        for(int j = 0; j < N; j++){
            vectres[i] += M[i][j] * vect[j];
        }
    }
}
// Calcula la norma del residuo ||Mx-b||
float ResidualNorm(float M[N][N], float x[N], float b[N]) {
    float res[N];
    for (int i = 0; i < N; i++) {
        res[i] = 0.0f;
        for (int j = 0; j < N; j++) {
            res[i] += M[i][j] * x[j];
        }
        res[i] -= b[i];
    }
    float sum = 0.0f;
    for (int i = 0; i < N; i++) {
        sum += res[i] * res[i];
    }
    return sqrt(sum);
}

/*
 hemos de dar error tras jacobi que da error. 
 Tras iter iteraciones obtenemos x_iter.
 Calculamos una iteración adicional x_next y el error = ||x_next - x_iter||_2 para que no seas valido.
*/
int Jacobi( float M[N][N], float b[N], float x_out[N], unsigned iter, float *error ){
    if(!DiagonalDom(M)) return 0;

    float x_curr[N]={0.0f};
    float x_new[N];

    // iter iteraciones normales
    for(unsigned k=0;k<iter;k++){
        for(int i=0;i<N;i++){
            float diag = M[i][i];
            if (fabsf(diag) < 1e-12f) return 0;
            float sum=0.0f;
            for(int j=0;j<N;j++)
                if(j!=i) sum += M[i][j]*x_curr[j];
            x_new[i] = (b[i]-sum)/diag;
        }
        for(int i=0;i<N;i++)
            x_curr[i]=x_new[i];
    }

    // Guardar solución tras 'iter' pasos
    for(int i=0;i<N;i++)
        x_out[i]=x_curr[i];

    // Una iteración extra para medir el “error” (diferencia próxima)
    float x_next[N];
    for(int i=0;i<N;i++){
        float diag = M[i][i];
        if (fabsf(diag) < 1e-12f) return 0;
        float sum=0.0f;
        for(int j=0;j<N;j++)
            if(j!=i) sum += M[i][j]*x_curr[j];
        x_next[i] = (b[i]-sum)/diag;
    }

    // Norma L2 del incremento próximo
    float diff = 0.0f;
    for(int i=0;i<N;i++){
        float d = x_next[i] - x_curr[i];
        diff += d*d;
    }
    if(error) *error = sqrtf(diff);

    return 1;
}

int main() {
    InitData();
    printf("V1 del 0 al 9 y del 256 al 265:\n");
    PrintVect(V1, 0, 10);
    PrintVect(V1, 256, 10);
    printf("V2 del 0 al 9 y del 256 al 265:\n");
    PrintVect(V2, 0, 10);
    PrintVect(V2, 256, 10);
    printf("V3 del 0 al 9 y del 256 al 265:\n");
    PrintVect(V3, 0, 10);
    PrintVect(V3, 256, 10);

    printf("Mat fila 0 y fila 100 del 0 al 9:\n");
    PrintRow(Mat, 0, 0, 10);
    PrintRow(Mat, 100, 0, 10);

    printf("MatDD fila 0 del 0 al 9 y fila 100 del 95 al 104:\n");
    PrintRow(MatDD, 0, 0, 10);
    PrintRow(MatDD, 100, 95, 10);

    printf("Infininorma de Mat = %.3f\n", Infininorm(Mat));
    printf("Norma de Frobenius de Mat = %.3f\n", NormFrobenius(Mat));
    printf("La matriz Mat %s diagonal dominante\n", DiagonalDom(Mat) ? "es" : "no es");

    printf("Infininorma de MatDD = %.3f\n", Infininorm(MatDD));
    printf("Norma de Frobenius de MatDD = %.3f\n", NormFrobenius(MatDD));
    printf("La matriz MatDD %s diagonal dominante\n", DiagonalDom(MatDD) ? "es" : "no es");

    printf("Escalar <V1,V2> = %f\n", Scalar(V1, V2));
    printf("Escalar <V1,V3> = %f\n", Scalar(V1, V3));
    printf("Escalar <V2,V3> = %f\n", Scalar(V2, V3));

    printf("Magnitud V1,V2 y V3 = %f %f %f\n", Magnitude(V1), Magnitude(V2), Magnitude(V3));

    printf("V1 y V2 %sson ortogonales\n", Ortogonal(V1, V2) ? "" : "no ");
    printf("V1 y V3 %sson ortogonales\n", Ortogonal(V1, V3) ? "" : "no ");
    printf("V2 y V3 %sson ortogonales\n", Ortogonal(V2, V3) ? "" : "no ");

    float V3x2[N];
    MultEscalar(V3, V3x2, 2.0f);
    printf("V3x2 del 0 al 9 y del 256 al 265:\n");
    PrintVect(V3x2, 0, 10);
    PrintVect(V3x2, 256, 10);

    float ProjV2onV3[N], ProjV1onV2[N];
    Projection(V2, V3, ProjV2onV3);
    Projection(V1, V2, ProjV1onV2);
    printf("Proyección de V2 sobre V3 (0-9):\n");
    PrintVect(ProjV2onV3, 0, 10);
    printf("Proyección de V1 sobre V2 (0-9):\n");
    PrintVect(ProjV1onV2, 0, 10);

    float MatV2[N];
    Matriu_x_Vector(Mat, V2, MatV2);
    printf("Mat*V2 (0-9):\n");
    PrintVect(MatV2, 0, 10);

    float sol1[N], sol1000[N];
    float err1=0.0f, err1000=0.0f;

    if (Jacobi(MatDD, V3, sol1, 1, &err1)){
        printf("Els elements 0 a 9 de la solució (1 iter) del sistema d'equacions són:\n");
        PrintVect(sol1,0,10);
        printf("Error %f\n", err1);
    } else {
        printf("La matriu MatDD no és diagonal dominant, no es pot aplicar Jacobi (1 iter)\n");
    }

    if (Jacobi(MatDD, V3, sol1000, 1000, &err1000)){
        printf("Els elements 0 a 9 de la solució (1000 iters) del sistema d'equacions són:\n");
        PrintVect(sol1000,0,10);
        printf("Error %f\n", err1000);
    } else {
        printf("La matriu MatDD no és diagonal dominant, no es pot aplicar Jacobi (1000 iters)\n");
    }

    if (!DiagonalDom(Mat)){
        printf("La matriu M no és diagonal dominant, no es pot aplicar Jacobi\n");
    } else {
        float solM[N], errM=0.0f;
        if (Jacobi(Mat, V3, solM, 1000, &errM)){
            PrintVect(solM,0,10);
            printf("Error %f\n", errM);
        } else {
            printf("La matriu M no és diagonal dominant, no es pot aplicar Jacobi\n");
        }
    }

    return 0;
}
