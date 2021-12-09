//Serial Code

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
FILE *ptr_file_1, *ptr_file_2;

#define RESET "\033[0m"
#define BOLDRED "\033[1m\033[31m"

#define PATH -1
#define NONE 0
#define UP 1
#define LEFT 2
#define DIAGONAL 3

int filelen1 = 0;
int filelen2 = 0;
int lenA, lenB;
char FASTA1[5000];			
char FASTA2[5000];	

void similarityScore(long long int i, long long int j, int *H, int *P, long long int *maxPos);
int matchMissmatchScore(long long int i, long long int j);
void backtrack(int *P, long long int maxPos);
void printMatrix(int *matrix);
void printPredecessorMatrix(int *matrix);
void generate(void);

long long int m = 11;
long long int n = 7;

int matchScore = 5;
int missmatchScore = -3;
int gapScore = -4;

char *a, *b;

int main(int argc, char *argv[])
{
    ptr_file_1 = fopen("s2.txt", "r");

    if (ptr_file_1 == NULL)
    {
        printf("Error opening 's2.txt'\n");
        system("PAUSE");
        exit(1);
    }

    ptr_file_2 = fopen("s1.txt", "r");

    if (ptr_file_2 == NULL)
    {
        printf("Error opening 's1.txt'\n");
        system("PAUSE");
        exit(1);
    }

    fgets(FASTA1, sizeof(FASTA1), ptr_file_1);
    fgets(FASTA2, sizeof(FASTA2), ptr_file_2);
    fclose(ptr_file_1);
    fclose(ptr_file_2);
	lenA = strlen(FASTA1);
	lenB = strlen(FASTA2);

    /*m = strtoll(argv[1], NULL, 10);
    n = strtoll(argv[2], NULL, 10);*/

    m = lenA;
    n = lenB;

#ifdef DEBUG
    printf("\nMatrix[%lld][%lld]\n", n, m);
#endif

    a = malloc(m * sizeof(char));
    b = malloc(n * sizeof(char));
    m++;
    n++;

    int *H;
    H = calloc(m * n, sizeof(int));
    int *P;
    P = calloc(m * n, sizeof(int));

    generate();

    long long int maxPos = 0;
    long long int i, j;
    double initialTime = omp_get_wtime();

    for (i = 1; i < n; i++)
    {
        for (j = 1; j < m; j++)
        {
            similarityScore(i, j, H, P, &maxPos);
        }
    }

    backtrack(P, maxPos);

    double finalTime = omp_get_wtime();
    printf("\nElapsed time: %f\n\n", finalTime - initialTime);

#ifdef DEBUG
    printf("\nSimilarity Matrix:\n");
    printMatrix(H);

    printf("\nPredecessor Matrix:\n");
    printPredecessorMatrix(P);
#endif

    free(H);
    free(P);
    free(a);
    free(b);

    return 0;
}

void similarityScore(long long int i, long long int j, int *H, int *P, long long int *maxPos)
{

    int up, left, diag;
    long long int index = m * i + j;
    up = H[index - m] + gapScore;
    left = H[index - 1] + gapScore;
    diag = H[index - m - 1] + matchMissmatchScore(i, j);

    int max = NONE;
    int pred = NONE;

    if (diag > max)
    {
        max = diag;
        pred = DIAGONAL;
    }

    if (up > max)
    {
        max = up;
        pred = UP;
    }

    if (left > max)
    {
        max = left;
        pred = LEFT;
    }

    H[index] = max;
    P[index] = pred;

    if (max > H[*maxPos])
    {
        *maxPos = index;
    }
}

int matchMissmatchScore(long long int i, long long int j)
{
    if (a[j - 1] == b[i - 1])
        return matchScore;
    else
        return missmatchScore;
}

void backtrack(int *P, long long int maxPos)
{
    long long int predPos;

    do
    {
        if (P[maxPos] == DIAGONAL)
            predPos = maxPos - m - 1;
        else if (P[maxPos] == UP)
            predPos = maxPos - m;
        else if (P[maxPos] == LEFT)
            predPos = maxPos - 1;
        P[maxPos] *= PATH;
        maxPos = predPos;
    } while (P[maxPos] != NONE);
}

void printMatrix(int *matrix)
{
    long long int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            printf("%d\t", matrix[m * i + j]);
        }
        printf("\n");
    }
}

void printPredecessorMatrix(int *matrix)
{
    long long int i, j, index;
    long long int countRed = 0;
    float similarityPercentage;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            index = m * i + j;
            if (matrix[index] < 0)
            {
                printf(BOLDRED);
                countRed += 1;
                if (matrix[index] == -UP)
                    printf("↑ ");
                else if (matrix[index] == -LEFT)
                    printf("← ");
                else if (matrix[index] == -DIAGONAL)
                    printf("↖ ");
                else
                    printf("- ");
                printf(RESET);
            }
            else
            {
                if (matrix[index] == UP)
                    printf("↑ ");
                else if (matrix[index] == LEFT)
                    printf("← ");
                else if (matrix[index] == DIAGONAL)
                    printf("↖ ");
                else
                    printf("- ");
            }
        }
        printf("\n");
    }
    similarityPercentage = (float) countRed / (i + j) * 100.0;
    printf("\n\n\t\tSimilarity Percentage = %.2f%%\n\n", similarityPercentage);
}

void generate()
{
    long long int i;
    for (i = 0; i < m; i++)
    {
        char aux = FASTA1[i];
        if (aux == 'a')
            a[i] = 'A';
        else if (aux == 'c')
            a[i] = 'C';
        else if (aux == 'g')
            a[i] = 'G';
        else if (aux == 't')
            a[i] = 'T';
    }

    for (i = 0; i < n; i++)
    {
        char aux = FASTA2[i];
        if (aux == 'a')
            b[i] = 'A';
        else if (aux == 'c')
            b[i] = 'C';
        else if (aux == 'g')
            b[i] = 'G';
        else if (aux == 't')
            b[i] = 'T';
    }
}

