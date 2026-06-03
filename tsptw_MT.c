/*
 * tsptw.c — Travelling Salesman Problem with Time Windows (TSPTW) Solver
 *
 * Maximizes visited cities (primary), minimizes tour length (secondary).
 *
 * Algorithm (per thread):
 *   1. Multi-start greedy construction (nearest / earliest-arrival / balanced)
 *   2. City insertion (greedy cheapest-feasible)
 *   3. 2-opt local search (limited segment, TW-aware)
 *   4. Or-opt / swap local search
 *   5. Perturbation (remove + reinsert) for ILS
 *   6. Tour rotation for best starting point
 *
 * Multi-threaded: NUM_THREADS independent solver threads share read-only
 * city data and merge results via a mutex-protected global best.
 *
 * Usage: tsptw.exe input.txt output.txt [time_limit_seconds]
 *
 * Compile: gcc -O3 -march=native -o tsptw tsptw.c -lm -lpthread
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

#ifdef _WIN32
#include <windows.h>
#endif

/* ===== Constants ===== */
#define MAX_N       55000
#define INF         0x7FFFFFFF
#define N_MODES     3
#define NUM_THREADS 4

/* ===== City data (read-only after init, shared across threads) ===== */
static int n;
static int cid[MAX_N], cx[MAX_N], cy[MAX_N], co[MAX_N], cc[MAX_N];

/* ===== Wall-clock Timing (thread-safe) ===== */
#ifdef _WIN32
static LARGE_INTEGER perf_freq, perf_start;
#else
static struct timespec ts_start;
#endif
static double tlim;

static double elapsed(void) {
#ifdef _WIN32
    LARGE_INTEGER now;
    QueryPerformanceCounter(&now);
    return (double)(now.QuadPart - perf_start.QuadPart) / (double)perf_freq.QuadPart;
#else
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - ts_start.tv_sec) + (now.tv_nsec - ts_start.tv_nsec) * 1e-9;
#endif
}

static int tup(void) { return elapsed() >= tlim; }

/* ===== Distance (Euclidean, rounded) ===== */
static inline int dist(int a, int b) {
    long long dx = cx[a] - cx[b], dy = cy[a] - cy[b];
    return (int)(sqrt((double)(dx * dx + dy * dy)) + 0.5);
}

static inline long long dist2(int a, int b) {
    long long dx = cx[a]-cx[b], dy = cy[a]-cy[b];
    return dx*dx + dy*dy;
}

/* ===== Spatial Grid (read-only after init, shared across threads) ===== */
static int gdim, cw, ch, gox, goy;
static int *gl, *gs, *gc_grid; /* list, start, count */

static void build_grid(void) {
    int mnx = cx[0], mxx = cx[0], mny = cy[0], mxy = cy[0];
    for (int i = 1; i < n; i++) {
        if (cx[i] < mnx) mnx = cx[i]; if (cx[i] > mxx) mxx = cx[i];
        if (cy[i] < mny) mny = cy[i]; if (cy[i] > mxy) mxy = cy[i];
    }
    gox = mnx; goy = mny;
    int rx = mxx - mnx + 1, ry = mxy - mny + 1;
    double area = (double)rx * ry;
    int cs = (int)(sqrt(area / (n / 25.0 + 1)) + 1);
    if (cs < 1) cs = 1;
    gdim = (rx > ry ? rx : ry) / cs + 1;
    if (gdim > 1000) gdim = 1000;
    if (gdim < 1) gdim = 1;
    cw = rx / gdim + 1; ch = ry / gdim + 1;
    if (cw < 1) cw = 1; if (ch < 1) ch = 1;

    int tc = gdim * gdim;
    gc_grid = (int *)calloc(tc, sizeof(int));
    gs = (int *)malloc(tc * sizeof(int));
    gl = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        int gx = (cx[i] - gox) / cw, gy = (cy[i] - goy) / ch;
        if (gx >= gdim) gx = gdim - 1; if (gy >= gdim) gy = gdim - 1;
        gc_grid[gy * gdim + gx]++;
    }
    gs[0] = 0;
    for (int i = 1; i < tc; i++) gs[i] = gs[i-1] + gc_grid[i-1];
    int *off = (int *)calloc(tc, sizeof(int));
    for (int i = 0; i < n; i++) {
        int gx = (cx[i] - gox) / cw, gy = (cy[i] - goy) / ch;
        if (gx >= gdim) gx = gdim - 1; if (gy >= gdim) gy = gdim - 1;
        int c = gy * gdim + gx;
        gl[gs[c] + off[c]++] = i;
    }
    free(off);
}

/* ===== Per-Thread Solver Context ===== */
typedef struct {
    /* Working tour */
    int *T;
    int K;
    int TL, TC;
    int *E;
    int *U;
    /* Thread-local best */
    int *BT;
    int BK, BL, BC;
    /* RNG state */
    unsigned long long rng_s;
    /* Scratch buffers (replacing static locals) */
    int *tmp;
    int *ne;
    int *mark;    /* for sel_starts */
    /* Thread identity */
    int thread_id;
    /* Starters assigned to this thread */
    int *starters;
    int num_starters;
} SolverCtx;

static SolverCtx *ctx_alloc(int tid, unsigned long long seed) {
    SolverCtx *ctx = (SolverCtx *)malloc(sizeof(SolverCtx));
    ctx->T    = (int *)malloc(MAX_N * sizeof(int));
    ctx->E    = (int *)malloc(MAX_N * sizeof(int));
    ctx->U    = (int *)calloc(MAX_N, sizeof(int));
    ctx->BT   = (int *)malloc(MAX_N * sizeof(int));
    ctx->tmp  = (int *)malloc(MAX_N * sizeof(int));
    ctx->ne   = (int *)malloc(MAX_N * sizeof(int));
    ctx->mark = (int *)calloc(MAX_N, sizeof(int));
    ctx->BK = 0; ctx->BL = INF; ctx->BC = INF;
    ctx->K = 0; ctx->TL = 0; ctx->TC = 0;
    ctx->rng_s = seed;
    ctx->thread_id = tid;
    ctx->starters = NULL;
    ctx->num_starters = 0;
    return ctx;
}

static void ctx_free(SolverCtx *ctx) {
    free(ctx->T);
    free(ctx->E);
    free(ctx->U);
    free(ctx->BT);
    free(ctx->tmp);
    free(ctx->ne);
    free(ctx->mark);
    if (ctx->starters) free(ctx->starters);
    free(ctx);
}

/* ===== Global Best (mutex-protected, shared across threads) ===== */
static int G_BT[MAX_N], G_BK = 0, G_BL = INF, G_BC = INF;
static pthread_mutex_t g_best_mutex = PTHREAD_MUTEX_INITIALIZER;

/* Merge thread-local best into global best. Returns 1 if global updated. */
static int merge_global_best(SolverCtx *ctx) {
    int updated = 0;
    pthread_mutex_lock(&g_best_mutex);
    if (ctx->BK > G_BK ||
        (ctx->BK == G_BK && ctx->BL < G_BL) ||
        (ctx->BK == G_BK && ctx->BL == G_BL && ctx->BC < G_BC)) {
        G_BK = ctx->BK; G_BL = ctx->BL; G_BC = ctx->BC;
        memcpy(G_BT, ctx->BT, ctx->BK * sizeof(int));
        updated = 1;
    }
    pthread_mutex_unlock(&g_best_mutex);
    return updated;
}

/* Load global best into thread-local best (for cross-pollination) */
static void load_global_best(SolverCtx *ctx) {
    pthread_mutex_lock(&g_best_mutex);
    if (G_BK > ctx->BK ||
        (G_BK == ctx->BK && G_BL < ctx->BL) ||
        (G_BK == ctx->BK && G_BL == ctx->BL && G_BC < ctx->BC)) {
        ctx->BK = G_BK; ctx->BL = G_BL; ctx->BC = G_BC;
        memcpy(ctx->BT, G_BT, G_BK * sizeof(int));
    }
    pthread_mutex_unlock(&g_best_mutex);
}

/* ===== RNG (LCG, per-thread) ===== */
static int rng(SolverCtx *ctx, int m) {
    ctx->rng_s = ctx->rng_s * 6364136223846793005ULL + 1442695040888963407ULL;
    return (int)((ctx->rng_s >> 33) % (unsigned)m);
}

/* ===== Evaluate tour: returns 1=feasible, 0=infeasible ===== */
static int eval(int *t, int k, int *len, int *comp, int *e) {
    if (k == 0) { *len = *comp = 0; return 1; }
    int ct = co[t[0]]; /* arrival at start = 0, effective = max(0, open) = open */
    if (ct > cc[t[0]]) return 0;
    e[0] = ct;
    int td = 0;
    for (int i = 1; i < k; i++) {
        int d = dist(t[i-1], t[i]);
        td += d; ct += d;
        int ev = ct > co[t[i]] ? ct : co[t[i]];
        if (ev > cc[t[i]]) return 0;
        e[i] = ev; ct = ev;
    }
    int d = dist(t[k-1], t[0]);
    td += d; ct += d;
    *len = td; *comp = ct;
    return 1;
}

/* Recompute arrival times from position 'from'. Returns 1=ok, 0=infeasible */
static int recomp(int *t, int k, int *e, int from, int safe_from) {
    int ct;
    if (from == 0) {
        ct = co[t[0]]; if (ct > cc[t[0]]) return 0;
        e[0] = ct; from = 1;
    } else {
        ct = e[from - 1];
    }
    for (int i = from; i < k; i++) {
        ct += dist(t[i-1], t[i]);
        int ev = ct > co[t[i]] ? ct : co[t[i]];
        if (ev > cc[t[i]]) return 0;
        if (i >= safe_from && ev == e[i]) return 1; /* safe early break */
        e[i] = ev; ct = ev;
    }
    return 1;
}

/* Update thread-local best if current is better. Returns 1 if updated. */
static int ubest(SolverCtx *ctx) {
    if (ctx->K > ctx->BK ||
        (ctx->K == ctx->BK && ctx->TL < ctx->BL) ||
        (ctx->K == ctx->BK && ctx->TL == ctx->BL && ctx->TC < ctx->BC)) {
        ctx->BK = ctx->K; ctx->BL = ctx->TL; ctx->BC = ctx->TC;
        memcpy(ctx->BT, ctx->T, ctx->K * sizeof(int));
        return 1;
    }
    return 0;
}

/* Load thread-local best solution back into working arrays */
static void lbest(SolverCtx *ctx) {
    memcpy(ctx->T, ctx->BT, ctx->BK * sizeof(int));
    ctx->K = ctx->BK; ctx->TL = ctx->BL; ctx->TC = ctx->BC;
    memset(ctx->U, 0, n * sizeof(int));
    for (int i = 0; i < ctx->K; i++) ctx->U[ctx->T[i]] = 1;
    eval(ctx->T, ctx->K, &ctx->TL, &ctx->TC, ctx->E);
}

/* ===== Starter Selection ===== */
static int sel_starts(SolverCtx *ctx, int *st) {
    int ns = 0;
    memset(ctx->mark, 0, n * sizeof(int));

    /* Limit per category based on instance size */
    int lim = n < 200 ? n : (n < 2000 ? 80 : (n < 10000 ? 40 : 25));

    /* 1) Earliest open times */
    for (int pass = 0; pass < lim; pass++) {
        int best = -1, bv = INF;
        for (int i = 0; i < n; i++)
            if (!ctx->mark[i] && co[i] < bv) { bv = co[i]; best = i; }
        if (best < 0) break;
        st[ns++] = best; ctx->mark[best] = 1;
    }
    /* 2) Tightest windows */
    for (int pass = 0; pass < lim; pass++) {
        int best = -1, bv = INF;
        for (int i = 0; i < n; i++) {
            int w = cc[i] - co[i];
            if (!ctx->mark[i] && w < bv) { bv = w; best = i; }
        }
        if (best < 0) break;
        st[ns++] = best; ctx->mark[best] = 1;
    }
    /* 3) Random sample */
    for (int pass = 0; pass < lim; pass++) {
        int c = rng(ctx, n);
        if (!ctx->mark[c]) { st[ns++] = c; ctx->mark[c] = 1; }
    }
    return ns;
}

/* ===== Greedy Construction =====
 * mode 0: nearest feasible (minimize distance)
 * mode 1: earliest effective arrival (minimize time consumption)
 * mode 2: balanced (distance + effective arrival)
 */
static void greedy(SolverCtx *ctx, int start, int mode) {
    memset(ctx->U, 0, n * sizeof(int));
    ctx->T[0] = start; ctx->K = 1; ctx->U[start] = 1;
    int ct = co[start]; ctx->E[0] = ct;

    int mcell = cw < ch ? cw : ch;
    if (mcell < 1) mcell = 1;

    while (ctx->K < n) {
        int last = ctx->T[ctx->K - 1];
        int gx0 = (cx[last] - gox) / cw, gy0 = (cy[last] - goy) / ch;
        if (gx0 >= gdim) gx0 = gdim - 1; if (gx0 < 0) gx0 = 0;
        if (gy0 >= gdim) gy0 = gdim - 1; if (gy0 < 0) gy0 = 0;

        int best = -1;
        long long bscore = 0x7FFFFFFFFFFFFFFFLL;
        for (int r = 0; r <= gdim; r++) {
            /* Early termination: lower bound on distance for ring r */
            if (r >= 2 && best >= 0) {
                int lb = (r - 1) * mcell;
                if (mode == 0 && (long long)lb * lb > bscore) break;
                if (mode == 1 && ct + lb > bscore) break;
                if (mode == 2 && lb + ct + lb > bscore) break;
            }
            for (int dy = -r; dy <= r; dy++) {
                int gy = gy0 + dy;
                if (gy < 0 || gy >= gdim) continue;
                for (int dx = -r; dx <= r; dx++) {
                    if (r > 0 && abs(dx) < r && abs(dy) < r) continue;
                    int gx = gx0 + dx;
                    if (gx < 0 || gx >= gdim) continue;
                    int cell = gy * gdim + gx;
                    for (int j = gs[cell]; j < gs[cell] + gc_grid[cell]; j++) {
                        int c = gl[j];
                        if (ctx->U[c]) continue;
                        
                        long long score;
                        if (mode == 0) {
                            score = dist2(last, c);
                            if (score >= bscore) continue;
                            int d = dist(last, c);
                            int arr = ct + d;
                            int ev = arr > co[c] ? arr : co[c];
                            if (ev > cc[c]) continue;
                        } else {
                            int d = dist(last, c);
                            int arr = ct + d;
                            int ev = arr > co[c] ? arr : co[c];
                            if (ev > cc[c]) continue;
                            if (mode == 1) score = ev;
                            else score = d + ev;
                        }
                        if (score < bscore) { bscore = score; best = c; }
                    }
                }
            }
        }
        if (best < 0) break;

        int d = dist(last, best);
        int arr = ct + d;
        int ev = arr > co[best] ? arr : co[best];
        ctx->T[ctx->K] = best; ctx->E[ctx->K] = ev; ctx->U[best] = 1; ct = ev; ctx->K++;
    }
    eval(ctx->T, ctx->K, &ctx->TL, &ctx->TC, ctx->E);
}

/* ===== City Insertion Pass =====
 * For each unvisited city, find cheapest feasible insertion position.
 * Returns number of cities inserted.
 */
static int ins_pass(SolverCtx *ctx) {
    int inserted = 0;
    for (int c = 0; c < n && !tup(); c++) {
        if (ctx->U[c]) continue;
        int bpos = -1;
        long long bcost = 0x7FFFFFFFFFFFFFFFLL;
        for (int j = 1; j <= ctx->K; j++) {
            int prev = ctx->T[j-1];
            int next = (j < ctx->K) ? ctx->T[j] : ctx->T[0];
            long long cost = dist2(prev, c) + dist2(c, next) - dist2(prev, next);
            if (cost >= bcost) continue;

            /* TW feasibility check */
            int d1 = dist(prev, c);
            int arr = ctx->E[j-1] + d1;
            int ev = arr > co[c] ? arr : co[c];
            if (ev > cc[c]) continue;

            int d2 = dist(c, next);
            if (j < ctx->K) {
                int arr2 = ev + d2;
                int ev2 = arr2 > co[next] ? arr2 : co[next];
                if (ev2 > cc[next]) continue;
                /* Propagate if effective time at next city increased */
                if (ev2 > ctx->E[j]) {
                    int ok = 1, ct2 = ev2;
                    for (int p = j + 1; p < ctx->K; p++) {
                        ct2 += dist(ctx->T[p-1], ctx->T[p]);
                        int ep = ct2 > co[ctx->T[p]] ? ct2 : co[ctx->T[p]];
                        if (ep > cc[ctx->T[p]]) { ok = 0; break; }
                        if (ep <= ctx->E[p]) break;
                        ct2 = ep;
                    }
                    if (!ok) continue;
                }
            }
            bcost = cost; bpos = j;
        }
        if (bpos >= 0) 
        {
            for (int i = ctx->K; i > bpos; i--) {
                ctx->T[i] = ctx->T[i-1];
                ctx->E[i] = ctx->E[i-1]; /* Fix: Shift E along with T */
            }
            ctx->T[bpos] = c; ctx->K++; ctx->U[c] = 1;
            /* Incremental update: safe to break after the inserted position */
            recomp(ctx->T, ctx->K, ctx->E, bpos, bpos + 1);
            inserted++;
        }
    }
    if (inserted > 0) {
        eval(ctx->T, ctx->K, &ctx->TL, &ctx->TC, ctx->E);
    }
    return inserted;
}

/* ===== 2-opt Pass =====
 * Try reversing segments [i..j] if it shortens the tour.
 * Limited segment length for large tours.
 */
static int two_opt_pass(SolverCtx *ctx) {
    int improved = 0;
    int maxseg = ctx->K < 300 ? ctx->K : 300;

    /* Randomize starting position so each thread explores different neighborhoods */
    int range = ctx->K - 2;
    if (range < 1) return 0;
    int offset = rng(ctx, range);

    for (int ii = 0; ii < range && !tup(); ii++) {
        int i = 1 + (ii + offset) % range;
        int jlim = i + maxseg; if (jlim >= ctx->K) jlim = ctx->K - 1;
        for (int j = i + 1; j <= jlim; j++) {
            int a = ctx->T[i-1], b = ctx->T[i], c2 = ctx->T[j];
            int d = (j + 1 < ctx->K) ? ctx->T[j+1] : ctx->T[0];
            long long delta = dist2(a, c2) + dist2(b, d) - dist2(a, b) - dist2(c2, d);
            if (delta >= 0) continue;

            /* Reverse segment [i..j] */
            for (int l = i, r = j; l < r; l++, r--) {
                int tmp = ctx->T[l]; ctx->T[l] = ctx->T[r]; ctx->T[r] = tmp;
            }
            /* Check TW feasibility from position i. Safe to break after j. */
            if (recomp(ctx->T, ctx->K, ctx->E, i, j + 1)) {
                int actual_delta = dist(a, c2) + dist(b, d) - dist(a, b) - dist(c2, d);
                ctx->TL += actual_delta;
                ctx->TC = ctx->E[ctx->K-1] + dist(ctx->T[ctx->K-1], ctx->T[0]);
                improved = 1;
            } else {
                /* Undo reversal */
                for (int l = i, r = j; l < r; l++, r--) {
                    int tmp = ctx->T[l]; ctx->T[l] = ctx->T[r]; ctx->T[r] = tmp;
                }
                eval(ctx->T, ctx->K, &ctx->TL, &ctx->TC, ctx->E); /* full restore */
            }
        }
    }
    return improved;
}

/* ===== Swap Pass =====
 * Try swapping adjacent cities for distance improvement.
 */
static int swap_pass(SolverCtx *ctx) {
    int improved = 0;
    int range = ctx->K - 2;
    if (range < 1) return 0;
    int offset = rng(ctx, range);

    for (int ii = 0; ii < range && !tup(); ii++) {
        int i = 1 + (ii + offset) % range;
        int a = ctx->T[i-1], b = ctx->T[i], c2 = ctx->T[i+1];
        int d = (i + 2 < ctx->K) ? ctx->T[i+2] : ctx->T[0];
        long long delta = dist2(a, c2) + dist2(b, d) - dist2(a, b) - dist2(c2, d);
        if (delta >= 0) continue;

        ctx->T[i] = c2; ctx->T[i+1] = b;
        if (recomp(ctx->T, ctx->K, ctx->E, i, i + 2)) {
            int actual_delta = dist(a, c2) + dist(b, d) - dist(a, b) - dist(c2, d);
            ctx->TL += actual_delta;
            ctx->TC = ctx->E[ctx->K-1] + dist(ctx->T[ctx->K-1], ctx->T[0]);
            improved = 1;
        } else {
            ctx->T[i] = b; ctx->T[i+1] = c2;
            eval(ctx->T, ctx->K, &ctx->TL, &ctx->TC, ctx->E);
        }
    }
    return improved;
}

/* ===== Or-opt Pass =====
 * Try relocating single cities to better positions (limited window).
 */
static int oropt_pass(SolverCtx *ctx) {
    if (ctx->K > 5000) return 0; /* too expensive for large tours */
    int improved = 0;
    for (int i = 1; i < ctx->K && !tup(); i++) {
        int c = ctx->T[i];
        int prev_i = ctx->T[i-1];
        int next_i = (i+1 < ctx->K) ? ctx->T[i+1] : ctx->T[0];
        long long rem_save = dist2(prev_i, c) + dist2(c, next_i) - dist2(prev_i, next_i);

        int bj = -1;
        long long bdelta = 0; /* must be strictly negative to accept */
        int window = 40;
        for (int jj = 1; jj <= ctx->K; jj++) {
            /* Skip positions that are same or adjacent to i */
            if (jj == i || jj == i + 1) continue;

            int pj = ctx->T[jj - 1];
            int nj = (jj < ctx->K) ? ctx->T[jj] : ctx->T[0];
            /* When computing insertion, skip if the position references the removed city */
            if (jj > i) {
                /* After removal, index shifts: T[jj] becomes T[jj-1] etc */
                if (jj == i + 1) continue;
                pj = ctx->T[jj - 1]; /* but T[i] will be removed... approximate */
                if (jj - 1 == i) pj = ctx->T[i - 1];
                nj = (jj < ctx->K) ? ctx->T[jj] : ctx->T[0];
            }

            long long ins_cost = dist2(pj, c) + dist2(c, nj) - dist2(pj, nj);
            long long delta = ins_cost - rem_save;
            if (delta < bdelta) { bdelta = delta; bj = jj; }

            window--;
            if (window <= 0 && bj >= 0) break;
        }
        if (bj >= 0) {
            /* Perform relocation by rebuilding the tour */
            int tk = 0;
            for (int p = 0; p < ctx->K; p++) {
                if (p == i) continue;
                int adj = (bj > i) ? bj - 1 : bj;
                if (tk == adj) ctx->tmp[tk++] = c;
                ctx->tmp[tk++] = ctx->T[p];
            }
            if (tk < ctx->K) ctx->tmp[tk++] = c; /* insert at end if needed */

            int nl, nc;
            if (eval(ctx->tmp, ctx->K, &nl, &nc, ctx->ne)) {
                if (nl < ctx->TL || (nl == ctx->TL && nc < ctx->TC)) {
                    memcpy(ctx->T, ctx->tmp, ctx->K * sizeof(int));
                    memcpy(ctx->E, ctx->ne, ctx->K * sizeof(int));
                    ctx->TL = nl; ctx->TC = nc;
                    improved = 1;
                }
            }
        }
    }
    return improved;
}

/* ===== Perturbation =====
 * Remove random cities and try reinserting.
 */
static void perturb(SolverCtx *ctx) {
    if (ctx->K < 5) return;
    /* Vary perturbation strength per thread via RNG */
    int base = ctx->K / 10;
    int extra = rng(ctx, ctx->K / 5 + 1);
    int nr = base + extra;
    if (nr < 3) nr = 3;
    if (nr > 120) nr = 120;

    for (int i = 0; i < nr && ctx->K > 2; i++) {
        int pos = 1 + rng(ctx, ctx->K - 1);
        ctx->U[ctx->T[pos]] = 0;
        for (int j = pos; j < ctx->K - 1; j++) ctx->T[j] = ctx->T[j+1];
        ctx->K--;
    }
    eval(ctx->T, ctx->K, &ctx->TL, &ctx->TC, ctx->E);
    ins_pass(ctx);
}

/* ===== Tour Rotation =====
 * Try all rotations to find the best feasible starting point.
 */
static void try_rotations_ctx(SolverCtx *ctx) {
    if (ctx->BK <= 1) return;
    int best_start = -1, best_comp = ctx->BC;

    for (int s = 0; s < ctx->BK; s++) {
        for (int i = 0; i < ctx->BK; i++) ctx->tmp[i] = ctx->BT[(s + i) % ctx->BK];
        int rl, rc;
        if (eval(ctx->tmp, ctx->BK, &rl, &rc, ctx->ne)) {
            /* Same length guaranteed for rotations, pick smallest completion */
            if (rc < best_comp) {
                best_comp = rc;
                best_start = s;
            }
        }
    }
    if (best_start > 0) {
        for (int i = 0; i < ctx->BK; i++) ctx->tmp[i] = ctx->BT[(best_start + i) % ctx->BK];
        memcpy(ctx->BT, ctx->tmp, ctx->BK * sizeof(int));
        ctx->BC = best_comp;
    }
}

/* ===== Thread Solver Entry Point ===== */
static void *solver_thread(void *arg) {
    SolverCtx *ctx = (SolverCtx *)arg;

    /* Phase 1: Multi-start greedy construction on assigned starters.
     * Vary time budget per thread so they desynchronize and start
     * Phase 2 at different times. */
    double p1_fracs[] = {0.30, 0.20, 0.15, 0.10};
    double p1_frac = (ctx->thread_id < 4) ? p1_fracs[ctx->thread_id] : 0.20;
    double clim = elapsed() + tlim * p1_frac;

    for (int s = 0; s < ctx->num_starters && elapsed() < clim; s++) {
        for (int m = 0; m < N_MODES; m++) {
            greedy(ctx, ctx->starters[s], m);
            ubest(ctx);
            if (elapsed() >= clim) break;
        }
    }

    /* Phase 2: Iterative improvement — fully independent per thread.
     * No cross-pollination: each thread explores its own neighborhood
     * starting from its own Phase 1 best. Results merge only at the end. */
    lbest(ctx);
    int stale = 0;
    int iter = 0;
    /* Each thread uses a different stale threshold for perturbation */
    int stale_limit = 2 + (ctx->thread_id % 3); /* 2, 3, 4, 2 */

    while (!tup()) {
        int pk = ctx->BK, pl = ctx->BL;

        ins_pass(ctx);
        ubest(ctx);

        if (!tup()) { two_opt_pass(ctx); ubest(ctx); }
        if (!tup()) { swap_pass(ctx); ubest(ctx); }
        if (!tup()) { oropt_pass(ctx); ubest(ctx); }
        if (!tup()) { ins_pass(ctx); ubest(ctx); }

        if (ctx->BK == pk && ctx->BL == pl) {
            stale++;
            if (stale >= stale_limit) {
                /* Try a fresh random start */
                int s = rng(ctx, n);
                greedy(ctx, s, rng(ctx, N_MODES));
                ins_pass(ctx);
                ubest(ctx);
                lbest(ctx);
                stale = 0;
            } else {
                lbest(ctx);
                perturb(ctx);
                ubest(ctx);
                lbest(ctx);
            }
        } else {
            stale = 0;
            lbest(ctx);
        }

        iter++;
    }

    /* Final merge: push thread-local best into global best */
    merge_global_best(ctx);
    fprintf(stderr, "  Thread %d finished: local_best=%d cities, len=%d, iters=%d\n",
            ctx->thread_id, ctx->BK, ctx->BL, iter);
    return NULL;
}

/* ===== I/O ===== */
static void read_input(const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f) { fprintf(stderr, "Cannot open %s\n", fn); exit(1); }
    n = 0;
    int id, x, y, o, c;
    while (fscanf(f, "%d %d %d %d %d", &id, &x, &y, &o, &c) == 5) {
        cid[n] = id; cx[n] = x; cy[n] = y; co[n] = o; cc[n] = c;
        n++;
        if (n >= MAX_N) break;
    }
    fclose(f);
}

static void write_output(const char *fn) {
    FILE *f = fopen(fn, "w");
    if (!f) { fprintf(stderr, "Cannot write %s\n", fn); exit(1); }
    fprintf(f, "%d %d %d\n", G_BK, G_BL, G_BC);
    for (int i = 0; i < G_BK; i++)
        fprintf(f, "%d\n", cid[G_BT[i]]);
    fprintf(f, "\n");
    fclose(f);
    fprintf(stderr, "Result: %d cities, length %d, completion %d (%.1fs, %d threads)\n",
            G_BK, G_BL, G_BC, elapsed(), NUM_THREADS);
}

/* ===== Entry Point ===== */
int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s input.txt output.txt [time_limit]\n", argv[0]);
        return 1;
    }

    /* Initialize wall-clock timer */
#ifdef _WIN32
    QueryPerformanceFrequency(&perf_freq);
    QueryPerformanceCounter(&perf_start);
#else
    clock_gettime(CLOCK_MONOTONIC, &ts_start);
#endif

    tlim = argc > 3 ? atof(argv[3]) : 55.0;
    unsigned long long base_seed = (unsigned long long)time(NULL) ^ 0xABCDEF0123456789ULL;

    read_input(argv[1]);
    fprintf(stderr, "Loaded %d cities. Time limit: %.0fs. Threads: %d\n", n, tlim, NUM_THREADS);

    if (n == 0) {
        FILE *f = fopen(argv[2], "w");
        fprintf(f, "0 0 0\n\n");
        fclose(f);
        return 0;
    }

    /* Build spatial grid (single-threaded, read-only after this) */
    build_grid();

    /* Generate all starters using a temporary context */
    SolverCtx *tmp_ctx = ctx_alloc(-1, base_seed);
    int *all_starters = (int *)malloc(MAX_N * sizeof(int));
    int total_starters = sel_starts(tmp_ctx, all_starters);
    ctx_free(tmp_ctx);

    /* Create threads and distribute starters round-robin */
    pthread_t threads[NUM_THREADS];
    SolverCtx *ctxs[NUM_THREADS];

    for (int i = 0; i < NUM_THREADS; i++) {
        unsigned long long seed = base_seed ^ ((unsigned long long)(i + 1) * 0x9E3779B97F4A7C15ULL);
        ctxs[i] = ctx_alloc(i, seed);

        /* Count starters for this thread */
        int cnt = 0;
        for (int s = i; s < total_starters; s += NUM_THREADS) cnt++;
        ctxs[i]->starters = (int *)malloc(cnt * sizeof(int));
        ctxs[i]->num_starters = cnt;
        int idx = 0;
        for (int s = i; s < total_starters; s += NUM_THREADS)
            ctxs[i]->starters[idx++] = all_starters[s];
    }
    free(all_starters);

    /* Launch all solver threads */
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_create(&threads[i], NULL, solver_thread, ctxs[i]);
    }

    /* Wait for all threads to finish */
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    /* Phase 3: Try all rotations for best start on global best */
    if (G_BK <= 20000) {
        SolverCtx *rot_ctx = ctx_alloc(-1, base_seed);
        rot_ctx->BK = G_BK; rot_ctx->BL = G_BL; rot_ctx->BC = G_BC;
        memcpy(rot_ctx->BT, G_BT, G_BK * sizeof(int));
        try_rotations_ctx(rot_ctx);
        /* Copy back */
        G_BK = rot_ctx->BK; G_BL = rot_ctx->BL; G_BC = rot_ctx->BC;
        memcpy(G_BT, rot_ctx->BT, G_BK * sizeof(int));
        ctx_free(rot_ctx);
    }

    /* Cleanup */
    for (int i = 0; i < NUM_THREADS; i++) ctx_free(ctxs[i]);

    write_output(argv[2]);

    /* Free grid */
    free(gl); free(gs); free(gc_grid);

    return 0;
}
