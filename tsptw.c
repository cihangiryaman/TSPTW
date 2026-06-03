/*
 * tsptw.c — Travelling Salesman Problem with Time Windows (TSPTW) Solver
 *
 * Maximizes visited cities (primary), minimizes tour length (secondary).
 *
 * Algorithm:
 *   1. Multi-start greedy construction (nearest / earliest-arrival / balanced)
 *   2. City insertion (greedy cheapest-feasible)
 *   3. 2-opt local search (limited segment, TW-aware)
 *   4. Or-opt / swap local search
 *   5. Perturbation (remove + reinsert) for ILS
 *   6. Tour rotation for best starting point
 *
 * Usage: tsptw.exe input.txt output.txt [time_limit_seconds]
 *
 * Compile: gcc -O3 -march=native -flto -o tsptw tsptw.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* ===== Constants ===== */
#define MAX_N   55000
#define INF     0x7FFFFFFF
#define N_MODES 3

/* ===== City data (struct of arrays for cache performance) ===== */
static int n;
static int cid[MAX_N], cx[MAX_N], cy[MAX_N], co[MAX_N], cc[MAX_N];

/* ===== RNG (LCG) ===== */
static unsigned long long rng_s;
static int rng(int m) {
    rng_s = rng_s * 6364136223846793005ULL + 1442695040888963407ULL;
    return (int)((rng_s >> 33) % (unsigned)m);
}

/* ===== Timing ===== */
static clock_t clk0;
static double tlim;
static double elapsed(void) { return (double)(clock() - clk0) / CLOCKS_PER_SEC; }
static int tup(void) { return elapsed() >= tlim; }

/* ===== Distance (Euclidean, rounded) ===== */
static inline int dist(int a, int b) {
    long long dx = cx[a] - cx[b], dy = cy[a] - cy[b];
    return (int)(sqrt((double)(dx * dx + dy * dy)) + 0.5);
}

/* ===== Spatial Grid ===== */
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

/* ===== Tour State ===== */
static int T[MAX_N], K;       /* working tour + count */
static int TL, TC;            /* tour length, completion time */
static int E[MAX_N];          /* effective arrival times */
static int U[MAX_N];          /* in-tour flags */

/* ===== Best Solution ===== */
static int BT[MAX_N], BK, BL, BC;

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
static int recomp(int *t, int k, int *e, int from) {
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
        if (ev == e[i]) return 1; /* no change from here on */
        e[i] = ev; ct = ev;
    }
    return 1;
}

/* Update best if current is better. Returns 1 if updated. */
static int ubest(void) {
    if (K > BK || (K == BK && TL < BL) || (K == BK && TL == BL && TC < BC)) {
        BK = K; BL = TL; BC = TC;
        memcpy(BT, T, K * sizeof(int));
        return 1;
    }
    return 0;
}

/* Load best solution back into working arrays */
static void lbest(void) {
    memcpy(T, BT, BK * sizeof(int));
    K = BK; TL = BL; TC = BC;
    memset(U, 0, n * sizeof(int));
    for (int i = 0; i < K; i++) U[T[i]] = 1;
    eval(T, K, &TL, &TC, E);
}

/* ===== Starter Selection ===== */
static int sel_starts(int *st) {
    int ns = 0;
    static int mark[MAX_N];
    memset(mark, 0, n * sizeof(int));

    /* Limit per category based on instance size */
    int lim = n < 200 ? n : (n < 2000 ? 80 : (n < 10000 ? 40 : 25));

    /* 1) Earliest open times */
    for (int pass = 0; pass < lim; pass++) {
        int best = -1, bv = INF;
        for (int i = 0; i < n; i++)
            if (!mark[i] && co[i] < bv) { bv = co[i]; best = i; }
        if (best < 0) break;
        st[ns++] = best; mark[best] = 1;
    }
    /* 2) Tightest windows */
    for (int pass = 0; pass < lim; pass++) {
        int best = -1, bv = INF;
        for (int i = 0; i < n; i++) {
            int w = cc[i] - co[i];
            if (!mark[i] && w < bv) { bv = w; best = i; }
        }
        if (best < 0) break;
        st[ns++] = best; mark[best] = 1;
    }
    /* 3) Random sample */
    for (int pass = 0; pass < lim; pass++) {
        int c = rng(n);
        if (!mark[c]) { st[ns++] = c; mark[c] = 1; }
    }
    return ns;
}

/* ===== Greedy Construction =====
 * mode 0: nearest feasible (minimize distance)
 * mode 1: earliest effective arrival (minimize time consumption)
 * mode 2: balanced (distance + effective arrival)
 */
static void greedy(int start, int mode) {
    memset(U, 0, n * sizeof(int));
    T[0] = start; K = 1; U[start] = 1;
    int ct = co[start]; E[0] = ct;

    int mcell = cw < ch ? cw : ch;
    if (mcell < 1) mcell = 1;

    while (K < n) {
        int last = T[K - 1];
        int gx0 = (cx[last] - gox) / cw, gy0 = (cy[last] - goy) / ch;
        if (gx0 >= gdim) gx0 = gdim - 1; if (gx0 < 0) gx0 = 0;
        if (gy0 >= gdim) gy0 = gdim - 1; if (gy0 < 0) gy0 = 0;

        int best = -1, bscore = INF;
        for (int r = 0; r <= gdim; r++) {
            /* Early termination: lower bound on distance for ring r */
            if (r >= 2 && best >= 0) {
                int lb = (r - 1) * mcell;
                if (mode == 0 && lb > bscore) break;
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
                        if (U[c]) continue;
                        int d = dist(last, c);
                        int arr = ct + d;
                        int ev = arr > co[c] ? arr : co[c];
                        if (ev > cc[c]) continue;
                        int score;
                        switch (mode) {
                            case 0: score = d; break;
                            case 1: score = ev; break;
                            default: score = d + ev; break;
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
        T[K] = best; E[K] = ev; U[best] = 1; ct = ev; K++;
    }
    eval(T, K, &TL, &TC, E);
}

/* ===== City Insertion Pass =====
 * For each unvisited city, find cheapest feasible insertion position.
 * Returns number of cities inserted.
 */
static int ins_pass(void) {
    int inserted = 0;
    for (int c = 0; c < n && !tup(); c++) {
        if (U[c]) continue;
        int bpos = -1, bcost = INF;
        for (int j = 1; j <= K; j++) {
            int prev = T[j-1];
            int next = (j < K) ? T[j] : T[0];
            int d1 = dist(prev, c), d2 = dist(c, next), d3 = dist(prev, next);
            int cost = d1 + d2 - d3;
            if (cost >= bcost) continue;

            /* TW feasibility check */
            int arr = E[j-1] + d1;
            int ev = arr > co[c] ? arr : co[c];
            if (ev > cc[c]) continue;

            if (j < K) {
                int arr2 = ev + d2;
                int ev2 = arr2 > co[next] ? arr2 : co[next];
                if (ev2 > cc[next]) continue;
                /* Propagate if effective time at next city increased */
                if (ev2 > E[j]) {
                    int ok = 1, ct2 = ev2;
                    for (int p = j + 1; p < K; p++) {
                        ct2 += dist(T[p-1], T[p]);
                        int ep = ct2 > co[T[p]] ? ct2 : co[T[p]];
                        if (ep > cc[T[p]]) { ok = 0; break; }
                        if (ep <= E[p]) break;
                        ct2 = ep;
                    }
                    if (!ok) continue;
                }
            }
            bcost = cost; bpos = j;
        }
        if (bpos >= 0) {
            for (int i = K; i > bpos; i--) T[i] = T[i-1];
            T[bpos] = c; K++; U[c] = 1;
            eval(T, K, &TL, &TC, E);
            inserted++;
        }
    }
    return inserted;
}

/* ===== 2-opt Pass =====
 * Try reversing segments [i..j] if it shortens the tour.
 * Limited segment length for large tours.
 */
static int two_opt_pass(void) {
    int improved = 0;
    int maxseg = K < 300 ? K : 300;

    for (int i = 1; i < K - 1 && !tup(); i++) {
        int jlim = i + maxseg; if (jlim >= K) jlim = K - 1;
        for (int j = i + 1; j <= jlim; j++) {
            int a = T[i-1], b = T[i], c2 = T[j];
            int d = (j + 1 < K) ? T[j+1] : T[0];
            int delta = dist(a, c2) + dist(b, d) - dist(a, b) - dist(c2, d);
            if (delta >= 0) continue;

            /* Reverse segment [i..j] */
            for (int l = i, r = j; l < r; l++, r--) {
                int tmp = T[l]; T[l] = T[r]; T[r] = tmp;
            }
            /* Check TW feasibility from position i */
            if (recomp(T, K, E, i)) {
                TL += delta;
                TC = E[K-1] + dist(T[K-1], T[0]);
                improved = 1;
            } else {
                /* Undo reversal */
                for (int l = i, r = j; l < r; l++, r--) {
                    int tmp = T[l]; T[l] = T[r]; T[r] = tmp;
                }
                eval(T, K, &TL, &TC, E); /* full restore */
            }
        }
    }
    return improved;
}

/* ===== Swap Pass =====
 * Try swapping adjacent cities for distance improvement.
 */
static int swap_pass(void) {
    int improved = 0;
    for (int i = 1; i < K - 1 && !tup(); i++) {
        int a = T[i-1], b = T[i], c2 = T[i+1];
        int d = (i + 2 < K) ? T[i+2] : T[0];
        int delta = dist(a, c2) + dist(b, d) - dist(a, b) - dist(c2, d);
        if (delta >= 0) continue;

        T[i] = c2; T[i+1] = b;
        if (recomp(T, K, E, i)) {
            TL += delta;
            TC = E[K-1] + dist(T[K-1], T[0]);
            improved = 1;
        } else {
            T[i] = b; T[i+1] = c2;
            eval(T, K, &TL, &TC, E);
        }
    }
    return improved;
}

/* ===== Or-opt Pass =====
 * Try relocating single cities to better positions (limited window).
 */
static int oropt_pass(void) {
    if (K > 5000) return 0; /* too expensive for large tours */
    int improved = 0;
    for (int i = 1; i < K && !tup(); i++) {
        int c = T[i];
        int prev_i = T[i-1];
        int next_i = (i+1 < K) ? T[i+1] : T[0];
        int rem_save = dist(prev_i, c) + dist(c, next_i) - dist(prev_i, next_i);

        int bj = -1, bdelta = 0; /* must be strictly negative to accept */
        int window = 40;
        for (int jj = 1; jj <= K; jj++) {
            /* Skip positions that are same or adjacent to i */
            if (jj == i || jj == i + 1) continue;

            int pj = T[jj - 1];
            int nj = (jj < K) ? T[jj] : T[0];
            /* When computing insertion, skip if the position references the removed city */
            if (jj > i) {
                /* After removal, index shifts: T[jj] becomes T[jj-1] etc */
                if (jj == i + 1) continue;
                pj = T[jj - 1]; /* but T[i] will be removed... approximate */
                if (jj - 1 == i) pj = T[i - 1];
                nj = (jj < K) ? T[jj] : T[0];
            }

            int ins_cost = dist(pj, c) + dist(c, nj) - dist(pj, nj);
            int delta = ins_cost - rem_save;
            if (delta < bdelta) { bdelta = delta; bj = jj; }

            window--;
            if (window <= 0 && bj >= 0) break;
        }
        if (bj >= 0) {
            /* Perform relocation by rebuilding the tour */
            static int tmp[MAX_N];
            int tk = 0;
            for (int p = 0; p < K; p++) {
                if (p == i) continue;
                int adj = (bj > i) ? bj - 1 : bj;
                if (tk == adj) tmp[tk++] = c;
                tmp[tk++] = T[p];
            }
            if (tk < K) tmp[tk++] = c; /* insert at end if needed */

            int nl, nc;
            static int ne[MAX_N];
            if (eval(tmp, K, &nl, &nc, ne)) {
                if (nl < TL || (nl == TL && nc < TC)) {
                    memcpy(T, tmp, K * sizeof(int));
                    memcpy(E, ne, K * sizeof(int));
                    TL = nl; TC = nc;
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
static void perturb(void) {
    if (K < 5) return;
    int nr = K / 8;
    if (nr < 3) nr = 3;
    if (nr > 80) nr = 80;

    for (int i = 0; i < nr && K > 2; i++) {
        int pos = 1 + rng(K - 1);
        U[T[pos]] = 0;
        for (int j = pos; j < K - 1; j++) T[j] = T[j+1];
        K--;
    }
    eval(T, K, &TL, &TC, E);
    ins_pass();
}

/* ===== Tour Rotation =====
 * Try all rotations to find the best feasible starting point.
 */
static void try_rotations(void) {
    if (BK <= 1) return;
    static int rt[MAX_N], re[MAX_N];
    int best_start = -1, best_comp = BC;

    for (int s = 0; s < BK; s++) {
        for (int i = 0; i < BK; i++) rt[i] = BT[(s + i) % BK];
        int rl, rc;
        if (eval(rt, BK, &rl, &rc, re)) {
            /* Same length guaranteed for rotations, pick smallest completion */
            if (rc < best_comp) {
                best_comp = rc;
                best_start = s;
            }
        }
    }
    if (best_start > 0) {
        static int tmp[MAX_N];
        for (int i = 0; i < BK; i++) tmp[i] = BT[(best_start + i) % BK];
        memcpy(BT, tmp, BK * sizeof(int));
        BC = best_comp;
    }
}

/* ===== Main Solver ===== */
static void solve(void) {
    BK = 0; BL = INF; BC = INF;
    build_grid();

    /* Phase 1: Multi-start greedy construction */
    static int starters[MAX_N];
    int ns = sel_starts(starters);
    double clim = elapsed() + tlim * 0.25;

    for (int s = 0; s < ns && elapsed() < clim; s++) {
        for (int m = 0; m < N_MODES; m++) {
            greedy(starters[s], m);
            ubest();
            if (elapsed() >= clim) break;
        }
    }

    /* Phase 2: Iterative improvement */
    lbest();
    int stale = 0;

    while (!tup()) {
        int pk = BK, pl = BL;

        ins_pass();
        ubest();

        if (!tup()) { two_opt_pass(); ubest(); }
        if (!tup()) { swap_pass(); ubest(); }
        if (!tup()) { oropt_pass(); ubest(); }
        if (!tup()) { ins_pass(); ubest(); }

        if (BK == pk && BL == pl) {
            stale++;
            if (stale >= 3) {
                /* Try a fresh start */
                int s = rng(n);
                greedy(s, rng(N_MODES));
                /* Run insertion on fresh greedy */
                ins_pass();
                ubest();
                lbest();
                stale = 0;
            } else {
                lbest();
                perturb();
                ubest();
                lbest();
            }
        } else {
            stale = 0;
            lbest();
        }
    }

    /* Phase 3: Try all rotations for best start */
    if (BK <= 20000) try_rotations();
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
    fprintf(f, "%d %d %d\n", BK, BL, BC);
    for (int i = 0; i < BK; i++)
        fprintf(f, "%d\n", cid[BT[i]]);
    fprintf(f, "\n");
    fclose(f);
    fprintf(stderr, "Result: %d cities, length %d, completion %d (%.1fs)\n",
            BK, BL, BC, elapsed());
}

/* ===== Entry Point ===== */
int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s input.txt output.txt [time_limit]\n", argv[0]);
        return 1;
    }
    clk0 = clock();
    tlim = argc > 3 ? atof(argv[3]) : 55.0;
    rng_s = (unsigned long long)time(NULL) ^ 0xABCDEF0123456789ULL;

    read_input(argv[1]);
    fprintf(stderr, "Loaded %d cities. Time limit: %.0fs\n", n, tlim);

    if (n == 0) {
        FILE *f = fopen(argv[2], "w");
        fprintf(f, "0 0 0\n\n");
        fclose(f);
        return 0;
    }

    solve();
    write_output(argv[2]);
    return 0;
}
