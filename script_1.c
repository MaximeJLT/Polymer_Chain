#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

typedef struct {
    int N;          // nombre de segments
    double b;       // longueur de segment
    int T;          // nombre de conformations
    int nb_bins;    // nombre de bins pour l'histogramme
    unsigned int seed;
} Params;

static void usage(const char *p) {
    fprintf(stderr,
        "Usage: %s [-N <int>] [-b <double>] [-T <int>] [-bins <int>] [-s <seed>]\n"
        "Defaults: N=100, b=3.0, T=1000, bins=50, seed=time\n", p);
}

static inline double urand(void) {
    // U[0,1)
    return (double)rand() / ((double)RAND_MAX + 1.0);
}

int main(int argc, char **argv) {
    Params par = { .N = 100, .b = 3.0, .T = 1000, .nb_bins = 50, .seed = (unsigned)time(NULL) };

    // Parse CLI
    for (int i=1; i<argc; ++i) {
        if (!strcmp(argv[i], "-N") && i+1<argc) par.N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-b") && i+1<argc) par.b = atof(argv[++i]);
        else if (!strcmp(argv[i], "-T") && i+1<argc) par.T = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-bins") && i+1<argc) par.nb_bins = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-s") && i+1<argc) par.seed = (unsigned)strtoul(argv[++i], NULL, 10);
        else { usage(argv[0]); return 1; }
    }
    if (par.N <= 0 || par.T <= 0 || par.b <= 0.0 || par.nb_bins <= 0) {
        usage(argv[0]); return 1;
    }

    srand(par.seed);

    // Fichiers de sortie
    FILE *fxyz = fopen("polymere.xyz", "w");
    if (!fxyz) { perror("polymere.xyz"); return 1; }

    FILE *fqcsv = NULL; // on écrira l'histogramme après calcul

    // Tableaux pour Q et Q^2
    double *Q = (double*)malloc((size_t)par.T * sizeof(double));
    if (!Q) { fprintf(stderr, "Allocation Q failed\n"); fclose(fxyz); return 1; }

    // Accumulateurs pour QA_squared
    double sum_Q2 = 0.0;

    // Boucle sur les conformations
    for (int t=0; t<par.T; ++t) {
        // R1 = (0,0,0)
        double Rx = 0.0, Ry = 0.0, Rz = 0.0;
        // Ecrire en-tête XYZ pour la conformation t
        fprintf(fxyz, "%d\n", par.N + 1);
        fprintf(fxyz, "conformations %d\n", t+1);
        // Monomère 1
        fprintf(fxyz, "C %8.3f %8.3f %8.3f\n", Rx, Ry, Rz);

        // Somme cumulative des vecteurs de liaison (comme sum_b dans MATLAB)
        double sx = 0.0, sy = 0.0, sz = 0.0;

        for (int i=2; i<=par.N+1; ++i) {
            // theta ~ U[0, 2pi), phi = acos(2u-1) (uniforme sur la sphère)
            double theta = 2.0 * M_PI * urand();
            double u = urand();
            double phi = acos(2.0*u - 1.0);

            // Coordonnées cartésiennes du bond vector (norme b)
            double x = par.b * (sin(phi) * cos(theta));
            double y = par.b * (sin(phi) * sin(theta));
            double z = par.b * cos(phi);

            sx += x; sy += y; sz += z;

            // R(i) = R(1) + sum_b  (R1 est (0,0,0))
            Rx = sx; Ry = sy; Rz = sz;

            // Ecrire le point
            fprintf(fxyz, "C %8.3f %8.3f %8.3f\n", Rx, Ry, Rz);
        }

        // Q = || R(N+1) - R(1) || = || sum_b ||
        double Qnorm = sqrt(sx*sx + sy*sy + sz*sz);
        Q[t] = Qnorm;
        sum_Q2 += Qnorm * Qnorm;
    }

    fclose(fxyz);

    // Moyenne quadratique de Q
    double QA_squared = sum_Q2 / (double)par.T;
    printf("La valeur calculee de la moyenne quadratique de Q est : %.6f\n", QA_squared);

    // Valeur theorique
    double QA_theorique = (double)par.N * par.b * par.b;
    printf("La valeur theorique de la moyenne quadratique de Q est : %.6f\n", QA_theorique);

    // Histogramme de Q comme histcounts(Q, edges) avec edges = linspace(min, max, nb_bins+1)
    // 1) min & max
    double qmin = Q[0], qmax = Q[0];
    for (int t=1; t<par.T; ++t) {
        if (Q[t] < qmin) qmin = Q[t];
        if (Q[t] > qmax) qmax = Q[t];
    }
    // éviter binWidth=0 si toutes les valeurs sont identiques (cas pathologique)
    if (qmax == qmin) qmax = qmin + 1e-12;

    int *counts = (int*)calloc((size_t)par.nb_bins, sizeof(int));
    if (!counts) { fprintf(stderr, "Allocation counts failed\n"); free(Q); return 1; }

    double binWidth = (qmax - qmin) / (double)par.nb_bins;

    for (int t=0; t<par.T; ++t) {
        int idx = (int)floor((Q[t] - qmin) / binWidth);
        if (idx < 0) idx = 0;
        if (idx >= par.nb_bins) idx = par.nb_bins - 1; // inclure la valeur max dans le dernier bin
        counts[idx]++;
    }

    // Probabilités P(Q) = counts / T, et bords (edges)
    fqcsv = fopen("hist_Q.csv", "w");
    if (!fqcsv) {
        perror("hist_Q.csv");
        free(counts); free(Q);
        return 1;
    }
    // En-tete
    fprintf(fqcsv, "edge_left,edge_right,count,probability\n");
    for (int k=0; k<par.nb_bins; ++k) {
        double left  = qmin + k * binWidth;
        double right = left + binWidth;
        double prob  = (double)counts[k] / (double)par.T;
        fprintf(fqcsv, "%.10f,%.10f,%d,%.10f\n", left, right, counts[k], prob);
    }
    fclose(fqcsv);

    printf("Fichiers ecrits : polymere.xyz (conformations), hist_Q.csv (histogramme de Q)\n");

    free(counts);
    free(Q);
    return 0;
}
