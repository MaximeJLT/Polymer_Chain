#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

typedef struct {
    double b;           // segment length
    int T;              // number of conformations
    int Nmin, Nmax;     // sweep range
    int Nstep;          // step
    unsigned int seed;  // RNG seed
} Params;

static void usage(const char* prog){
    fprintf(stderr,
        "Usage: %s [-b <double>] [-T <int>] [-Nmin <int>] [-Nmax <int>] [-Nstep <int>] [-s <seed>]\n"
        "Defaults: b=3.0, T=1000, Nmin=100, Nmax=1000, Nstep=100, seed=time\n", prog);
}

static inline double urand01(void){
    // U[0,1)
    return (double)rand() / ((double)RAND_MAX + 1.0);
}

// Least squares linear regression y ~ m x + c, with R^2
static void linreg(const double *x, const double *y, int n, double *m, double *c, double *R2){
    double sx=0.0, sy=0.0, sxx=0.0, sxy=0.0;
    for(int i=0;i<n;++i){
        sx  += x[i];
        sy  += y[i];
        sxx += x[i]*x[i];
        sxy += x[i]*y[i];
    }
    double nx = (double)n;
    double denom = (nx * sxx - sx * sx);
    if (fabs(denom) < 1e-15) {
        *m = 0.0; *c = sy/nx; *R2 = 0.0; return;
    }
    *m = (nx * sxy - sx * sy) / denom;
    *c = (sy - (*m) * sx) / nx;

    // R^2 = 1 - SSE/SST
    double ymean = sy / nx;
    double SSE = 0.0, SST = 0.0;
    for(int i=0;i<n;++i){
        double yi_hat = (*m) * x[i] + (*c);
        double err = y[i] - yi_hat;
        SSE += err * err;
        double dev = y[i] - ymean;
        SST += dev * dev;
    }
    *R2 = (SST > 0.0) ? (1.0 - SSE / SST) : 0.0;
}

int main(int argc, char **argv){
    Params P = {.b=3.0, .T=1000, .Nmin=100, .Nmax=1000, .Nstep=100, .seed=(unsigned)time(NULL)};

    // Parse CLI
    for(int i=1;i<argc;i++){
        if (!strcmp(argv[i], "-b") && i+1<argc)       P.b = atof(argv[++i]);
        else if (!strcmp(argv[i], "-T") && i+1<argc)  P.T = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-Nmin") && i+1<argc) P.Nmin = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-Nmax") && i+1<argc) P.Nmax = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-Nstep") && i+1<argc) P.Nstep = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-s") && i+1<argc)    P.seed = (unsigned)strtoul(argv[++i], NULL, 10);
        else { usage(argv[0]); return 1; }
    }
    if (P.T <= 0 || P.b <= 0.0 || P.Nmin <= 0 || P.Nmax < P.Nmin || P.Nstep <= 0){
        usage(argv[0]); return 1;
    }

    srand(P.seed);

    // Determine number of N points
    int K = (P.Nmax - P.Nmin) / P.Nstep + 1;
    double *Nvals = (double*)malloc((size_t)K * sizeof(double));
    double *QA2   = (double*)malloc((size_t)K * sizeof(double));
    double *QA2th = (double*)malloc((size_t)K * sizeof(double));
    if(!Nvals || !QA2 || !QA2th){
        fprintf(stderr, "Allocation failed.\n");
        free(Nvals); free(QA2); free(QA2th);
        return 1;
    }

    // Sweep over N = Nmin:Nstep:Nmax
    for(int k=0;k<K;++k){
        int N = P.Nmin + k * P.Nstep;
        Nvals[k] = (double)N;

        // Accumulate mean of Q^2 over T conformations
        double sum_Q2 = 0.0;

        for(int t=0; t<P.T; ++t){
            // sum of bond vectors for this conformation
            double sx=0.0, sy=0.0, sz=0.0;

            for(int i=0; i<N; ++i){
                // theta uniform in [0, 2pi), phi = acos(2u-1)
                double theta = 2.0 * M_PI * urand01();
                double u = urand01();
                double phi = acos(2.0 * u - 1.0);

                // bond vector components of length b
                double x = P.b * (sin(phi) * cos(theta));
                double y = P.b * (sin(phi) * sin(theta));
                double z = P.b * cos(phi);

                sx += x; sy += y; sz += z;
            }

            double Qnorm2 = sx*sx + sy*sy + sz*sz;
            sum_Q2 += Qnorm2;
        }

        QA2[k]   = sum_Q2 / (double)P.T;       // mean(Q.^2)
        QA2th[k] = (double)N * P.b * P.b;      // theory: N*b^2
    }

    // Linear regression QA2 vs N
    double slope=0.0, intercept=0.0, R2=0.0;
    linreg(Nvals, QA2, K, &slope, &intercept, &R2);
    double b_est = (slope > 0.0) ? sqrt(slope) : 0.0;

    // Print summary
    printf("=== FJC sweep results ===\n");
    printf("b (input) = %.6f, T = %d, N in [%d:%d:%d]\n", P.b, P.T, P.Nmin, P.Nstep, P.Nmax);
    printf("Linear fit QA_squared ~ slope*N + intercept\n");
    printf("slope = %.6f, intercept = %.6f, R^2 = %.6f\n", slope, intercept, R2);
    printf("Estimated segment length b = sqrt(slope) = %.6f\n", b_est);

    // Write CSV for plotting
    FILE *f = fopen("qa_vs_N.csv", "w");
    if(!f){
        perror("qa_vs_N.csv");
    } else {
        fprintf(f, "N,QA_squared,QA_theory\n");
        for(int k=0;k<K;++k){
            fprintf(f, "%.0f,%.10f,%.10f\n", Nvals[k], QA2[k], QA2th[k]);
        }
        fclose(f);
        printf("CSV written: qa_vs_N.csv (columns: N, QA_squared, QA_theory)\n");
    }

    free(Nvals); free(QA2); free(QA2th);
    return 0;
}
