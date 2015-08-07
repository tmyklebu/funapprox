/*
Tor Myklebust's function approximation heuristic
Copyright (C) 2014-2015 Tor Myklebust (tmyklebu@csclub.uwaterloo.ca)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <queue>
#include <algorithm>
#include <sys/types.h>
#include <sys/resource.h>
#include <gmpxx.h>
#include <vector>
extern "C" {
#include <QSopt_ex.h>
}
using namespace std;

#define FOR(i,n) for (int i=0;i<(signed int)(n);i++)

long long get_cpu_usecs() {
 rusage r;
 getrusage(RUSAGE_SELF, &r);
 return            r.ru_utime.tv_usec
     + 1000000LL * r.ru_utime.tv_sec;
}

long long program_start;

// Wrapper class for a QSopt_ex rational linear program.
struct lp_t {
  mpq_QSprob prob;
  int nvar;
  int nextrow;

  lp_t(int nvar, const float *lo, const float *up) : nvar(nvar) {
    nextrow = 0;
    prob = mpq_QScreate_prob("funapprox", QS_MIN);
    FOR(i,nvar) {
      mpq_class zero = 0;
      mpq_class l = lo[i], h = up[i];
      mpq_QSnew_col(prob, zero.get_mpq_t(), l.get_mpq_t(), h.get_mpq_t(), 0);
    }
  }

  ~lp_t() {
    mpq_QSfree_prob(prob);
  }

  mpq_class solve() {
    int stat;
    if (QSexact_solver(prob, 0, 0, 0, DUAL_SIMPLEX, &stat)) {
      printf("solve failed\n");
      throw ":(";
    }
    if (stat != QS_LP_OPTIMAL) {
      throw ":(";
    }
    mpq_t foo;
    mpq_init(foo);
    mpq_QSget_objval(prob, &foo);
    mpq_class ans = mpq_class(foo);
    mpq_clear(foo);
    return ans;
  }

  void set_objective(vector<mpq_class> obj) {
    FOR(i, nvar) mpq_QSchange_objcoef(prob, i, obj[i].get_mpq_t());
  }

  void introduce_row(vector<mpq_class> lhs, char kind, mpq_class rhs) {
    mpq_QSnew_row(prob, rhs.get_mpq_t(), kind, 0);
    FOR(i, nvar)
      mpq_QSchange_coef(prob, nextrow, i, lhs[i].get_mpq_t());
    nextrow++;
  }

  vector<mpq_class> get_current_solution() {
    mpq_t foo[nvar];
    FOR(i,nvar) mpq_init(foo[i]);
    mpq_QSget_x_array(prob, foo);
    vector<mpq_class> ans(nvar);
    FOR(i,nvar) ans[i] = mpq_class(foo[i]);
    FOR(i,nvar) mpq_clear(foo[i]);
    return ans;
  }
};

// Find all `x` such that `lb <= fma(x, y, z) < ub`.  The output is the
// interval `[lo, hi)`.
//
// The implementation of this function is quite savage.  It could profitably be
// replaced by a couple of binary searches.  However, this implementation works
// and I'm not too keen on touching it.
void reverse_fmaf(float &lo, float &hi, float y, float z, float lb, float ub) {
  if (y < 0) {
    reverse_fmaf(lo, hi, -y, -z, -nextafterf(ub, -1.0/0.0), -nextafterf(lb, -1.0/0.0));
    assert(fmaf(lo, y, z) < ub);
    assert(fmaf(hi, y, z) < lb);
    assert(fmaf(nextafterf(lo, -1.0/0.0), y, z) >= ub);
    assert(fmaf(nextafterf(hi, -1.0/0.0), y, z) >= lb);
    return;
  }

  lo = -FLT_MAX; hi = FLT_MAX;
  while (nextafterf(lo, hi) < hi) {
    float c = lo/2 + hi/2;   // This BS is bad for floats.
    float cc = fmaf(c, y, z);
    if (cc < lb) lo = nextafterf(c, 1.0/0.0);
    else hi = c;
  }
  lo = nextafterf(lo, -1.0/0.0);
  lo = nextafterf(lo, -1.0/0.0);
  assert(fmaf(lo, y, z) < lb);
  while (fmaf(lo, y, z) < lb) lo = nextafterf(lo, 1.0/0.0);
  hi = lo;
  #define STEPIT(step) while (hi + (float)step != hi && fmaf(hi+(float)step, y, z) < ub) hi = hi + (float)step;
  for (int zzz = 0; zzz < 3; zzz++) {
    STEPIT(0x1.0p+120)
    STEPIT(0x1.0p+110)
    STEPIT(0x1.0p+100)
    STEPIT(0x1.0p+90)
    STEPIT(0x1.0p+80)
    STEPIT(0x1.0p+70)
    STEPIT(0x1.0p+60)
    STEPIT(0x1.0p+50)
    STEPIT(0x1.0p+45)
    STEPIT(0x1.0p+40)
    STEPIT(0x1.0p+35)
    STEPIT(0x1.0p+30)
    STEPIT(0x1.0p+25)
    STEPIT(0x1.0p+20)
    STEPIT(0x1.0p+15)
    STEPIT(0x1.0p+10)
    STEPIT(0x1.0p+5)
    STEPIT(0x1.0p+0)
    STEPIT(0x1.0p-5)
    STEPIT(0x1.0p-10)
    STEPIT(0x1.0p-15)
    STEPIT(0x1.0p-20)
    STEPIT(0x1.0p-25)
    STEPIT(0x1.0p-30)
    STEPIT(0x1.0p-35)
    STEPIT(0x1.0p-40)
    STEPIT(0x1.0p-45)
    STEPIT(0x1.0p-50)
    STEPIT(0x1.0p-55)
    STEPIT(0x1.0p-60)
    STEPIT(0x1.0p-65)
    STEPIT(0x1.0p-70)
    STEPIT(0x1.0p-80)
    STEPIT(0x1.0p-90)
    STEPIT(0x1.0p-100)
    STEPIT(0x1.0p-110)
    STEPIT(0x1.0p-120)
  }
  #undef STEPIT
  while (fmaf(hi, y, z) < ub) hi = nextafterf(hi, 1.0/0.0);
  assert(fmaf(lo, y, z) >= lb);
  assert(fmaf(hi, y, z) >= ub);
  assert(fmaf(nextafterf(lo, -1.0/0.0), y, z) < lb);
  assert(fmaf(nextafterf(hi, -1.0/0.0), y, z) < ub);
}

// A node in a straight-line program that only has constants, variables, and
// fused multiply-adds.
struct expression {
  int tag;
  union {
    float val;
    int varno;
    struct { expression *a, *b, *c; };
  };

  void getbounds(const float *xl, const float *xu, float &lower, float &upper) {
    switch (tag) {
      case 0: lower = upper = val; return;
      case 1: lower = xl[varno]; upper = xu[varno]; return;
      case 2: {
        float al, au, bl, bu, cl, cu;
        a->getbounds(xl, xu, al, au);
        b->getbounds(xl, xu, bl, bu);
        c->getbounds(xl, xu, cl, cu);
        float v1 = fmaf(al, bl, cl);
        float v2 = fmaf(au, bl, cl);
        float v3 = fmaf(al, bu, cl);
        float v4 = fmaf(au, bu, cl);
        float v5 = fmaf(al, bl, cu);
        float v6 = fmaf(au, bl, cu);
        float v7 = fmaf(al, bu, cu);
        float v8 = fmaf(au, bu, cu);
        lower = upper = v1;
        lower = min(lower, v2); upper = max(upper, v2);
        lower = min(lower, v3); upper = max(upper, v3);
        lower = min(lower, v4); upper = max(upper, v4);
        lower = min(lower, v5); upper = max(upper, v5);
        lower = min(lower, v6); upper = max(upper, v6);
        lower = min(lower, v7); upper = max(upper, v7);
        lower = min(lower, v8); upper = max(upper, v8);
      } return;
      default: abort();
    }
  }

  static expression *con(float f) {
    expression *e = new expression();
    e->tag = 0;
    e->val = f;
    return e;
  }

  static expression *var(int k) {
    expression *e = new expression();
    e->tag = 1;
    e->varno = k;
    return e;
  }

  static expression *fma(expression *a, expression *b, expression *c)  {
    expression *e = new expression();
    e->tag = 2;
    e->a = a; e->b = b; e->c = c;
    return e;
  }
};

expression *constant_fold(expression *e) {
  switch (e->tag) {
    case 0: return new expression(*e);
    case 1: return new expression(*e);
    case 2: {
      expression *a2 = constant_fold(e->a), *b2 = constant_fold(e->b),
                 *c2 = constant_fold(e->c);
      if (a2->tag == 0 && b2->tag == 0 && c2->tag == 0)
        return expression::con(fmaf(a2->val, b2->val, c2->val));
      else return expression::fma(a2, b2, c2);
    }
    default: abort();
  }
}

expression *subs(expression *e, int varno, float value) {
  switch (e->tag) {
    case 0: return new expression(*e);
    case 1: {
      if (e->varno == varno) return expression::con(value);
      else return new expression(*e);
    } abort();
    case 2: return expression::fma(subs(e->a, varno, value),
                                   subs(e->b, varno, value),
                                   subs(e->c, varno, value));
    default: abort();
  }
}

float eval(expression *e, float *x) {
  switch (e->tag) {
    case 0: return e->val;
    case 1: return x[e->varno];
    case 2: return fmaf(eval(e->a, x), eval(e->b, x), eval(e->c, x));
    default: abort();
  }
}

void pretty_print(expression *e) {
  switch (e->tag) {
    case 0: printf("%a", e->val); return;
    case 1: printf("x%i", e->varno); return;
    case 2: printf("fma(");
       pretty_print(e->a); printf(", ");
       pretty_print(e->b); printf(", ");
       pretty_print(e->c); printf(")");
    return;
    default: abort();
  }
}




float half_ulp(float x) {
  x = fabs(x);
  return (nextafterf(x, 1.0/0.0) - x) / 2;
}

// Compute the range `[lower, upper)` of acceptable function values at `x`.
//
// Hack this function, and `main()` below, if you'd like to do something
// other than approximate `tan(x)` on `(10^{-4}, pi/4)` within `0.999` ulp.
static double ulps_wrong = 0.999;
void get_bounds(float x, float &lower, float &upper) {
  double d = tan((double)x);
  float ulp = 2 * half_ulp((float)d);
  double lo = d - ulps_wrong * ulp, up = d + ulps_wrong * ulp;
  lower = lo; upper = up;
  if (lower < lo) lower = nextafterf(lower, 1.0/0.0);
  if (upper > up) upper = nextafterf(upper, -1.0/0.0);
  upper = nextafterf(upper, 1.0/0.0);
}

// Given a straight-line program `e` and bounds
// `[clb_0, cub_0] x ... x [clb_{nvar-1}, cub_{nvar-1}]` on the unspecified
// coefficients in `e`, this function computes an upper bound on the roundoff
// error incurred when evaluating `e` in `float` arithmetic.
mpq_class get_max_roundoff(expression *e, const float *clb, const float *cub) {
  switch (e->tag) {
    case 0: return 0;
    case 1: return 0;
    case 2: {
      float alo, ahi, blo, bhi, clo, chi;
      float elo, ehi;
      e->a->getbounds(clb, cub, alo, ahi);
      e->b->getbounds(clb, cub, blo, bhi);
      e->c->getbounds(clb, cub, clo, chi);
      e->getbounds(clb, cub, elo, ehi);
      float abound = max(fabs(alo), fabs(ahi));
      float bbound = max(fabs(blo), fabs(bhi));
      float ebound = max(fabs(elo), fabs(ehi));
      mpq_class around = get_max_roundoff(e->a, clb, cub);
      mpq_class bround = get_max_roundoff(e->b, clb, cub);
      mpq_class cround = get_max_roundoff(e->c, clb, cub);
      return half_ulp(ebound) + cround + bbound * around + abound * bround;
    }
    default: abort();
  }
}

// Returns the (rational) coefficients of the linear polynomial represented by
// the straight-line program `e`; `poly` gets the degree-1 terms and `con` gets
// the constant term.
//
// Each `fma` encountered by this function must have the form
// `fma(constant, x, y)` or `fma(x, constant, y)` where `x` and `y` are
// arbitrary.  If an `fma` encountered does not have this form, this function
// calls `abort()`.  Further, if `e` contains a variable for which `poly` does
// not have a corresponding index, undefined behaviour results.
void get_linear_poly(expression *e, vector<mpq_class> &poly, mpq_class &con) {
  FOR(i, poly.size()) poly[i] = 0;
  con = 0;
  switch (e->tag) {
    case 0: con = e->val; return;
    case 1: poly[e->varno] = 1; return;
    case 2: {
      vector<mpq_class> cp(poly.size());
      mpq_class cc;
      get_linear_poly(e->c, cp, cc);
      float constant;
      if (e->a->tag == 0) {
        constant = e->a->val;
        get_linear_poly(e->b, poly, con);
      } else if (e->b->tag == 0) {
        constant = e->b->val;
        get_linear_poly(e->a, poly, con);
      } else abort();
      FOR(i,poly.size()) poly[i] = constant * poly[i] + cp[i];
      con = constant * con + cc;
      return;
    }
    default: abort();
  }
}

// Peel back fma(const, f, const) -> f using reverse_fma.
// Now we have an expression tree and some bounds on the value it must take.
// Convert the expression tree to an inequality pair.
void peel_bounded_expression(expression *&e, float &lower, float &upper) {
  while (e->tag == 2 && e->c->tag == 0
         && (e->a->tag == 0 || e->b->tag == 0)) {
    if (e->a->tag == 0 && e->b->tag == 0) abort();
    float lo2, up2;
    float cee = e->c->val, bee;
    if (e->a->tag == 0) bee = e->a->val;
    else if (e->b->tag == 0) bee = e->b->val;
    else abort();
    reverse_fmaf(lo2, up2, bee, cee, lower, upper);
    lower = lo2; upper = up2;
    if (e->a->tag == 0) e = e->b; else e = e->a;
  }
}

// This function does the dirty work of `gen_inequalities` below.
void gen_inequalities_inner(expression *e, float lower, float upper, int nvar,
    const float *clb, const float *cub, lp_t *lp) {
  peel_bounded_expression(e, lower, upper);
  upper = nextafterf(upper, -1.0/0.0);

  vector<mpq_class> poly(nvar); mpq_class con;
  get_linear_poly(e, poly, con);

  mpq_class err = get_max_roundoff(e, clb, cub);
  mpq_class rhslo = lower - con - err, rhsup = upper - con + err;

  lp->introduce_row(poly, 'G', rhslo);
  lp->introduce_row(poly, 'L', rhsup);
}

// Given a straight-line program `e`, the number of unspecified coefficients
// `nvar`, bounds `[clb_0, cub_0] x ... x [clb_{nvar-1}, cub_{nvar-1}]` on the
// unspecified coefficients, and a test point `x`, add two inequalities to the
// linear program `lp` that provably must be satisfied by the unspecified
// coefficients.  These inequalities are derived by considering evaluation at
// the test point `x`.
void gen_inequalities(expression *e, int nvar,
    const float *clb, const float *cub, float x, lp_t *lp) {
  expression *ee = constant_fold(subs(e, -1, x));
  float lower, upper;
  get_bounds(x, lower, upper);
  gen_inequalities_inner(ee, lower, upper, nvar, clb, cub, lp);
}

// Compute lower and upper bounds on variable `varno` implied by the polyhedron
// defined by the constraints of `lp`.  Outputs a range `[lower, upper)`.
void get_var_bounds(lp_t *lp, int varno, float &lower, float &upper) {
  vector<mpq_class> obj(lp->nvar);
  obj[varno] = 1; lp->set_objective(obj);
  mpq_class lo = lp->solve();
  obj[varno] = -1; lp->set_objective(obj);
  mpq_class up = -lp->solve();
  lower = lo.get_d(); // XXX Possible double-rounding anomaly?
  upper = up.get_d(); // XXX Possible double-rounding anomaly?
  while (lower < lo) lower = nextafterf(lower, 1.0/0.0);
  while (upper > up) upper = nextafterf(upper, -1.0/0.0);
  upper = nextafterf(upper, 1.0/0.0);
}

// Bias guesses toward the centre of the interval.
// There is probably a better distribution.
double randpt() {
  return (drand48() + drand48()) / 2;
}

// This is the body of the diving heuristic.
vector<float> dive(int nvar, const float *clb, const float *cub,
    const vector<pair<expression *, pair<float, float> > > &ineqs,
    const vector<float> &preferred, int tries = 8) {
  if (nvar == 0) return vector<float>();
  lp_t lp(nvar, clb, cub);
  FOR(i, ineqs.size()) {
    gen_inequalities_inner(ineqs[i].first, ineqs[i].second.first,
        ineqs[i].second.second, nvar, clb, cub, &lp);
  }
  float alo, ahi;
  get_var_bounds(&lp, nvar-1, alo, ahi);
  ahi = nextafterf(ahi, -1.0/0.0);

  double num_ulps = (cub[nvar-1] - clb[nvar-1])
      / half_ulp(max(fabs(clb[nvar-1]), fabs(cub[nvar-1])));
  printf("%a\n", num_ulps);

  if (num_ulps > tries) {
    FOR(zzz, tries) {
      float f = randpt() * (ahi - alo) + alo;
      if (zzz == 0 && preferred[nvar-1] >= alo && preferred[nvar-1] <= ahi)
        f = preferred[nvar-1];
      vector<pair<expression *, pair<float, float> > > new_ineqs;
      FOR(i, ineqs.size()) new_ineqs.push_back(make_pair(
          constant_fold(subs(ineqs[i].first, nvar-1, f)),
          ineqs[i].second));
      try {
        vector<float> ans = dive(nvar-1, clb, cub, new_ineqs, preferred);
        ans.push_back(f);
        return ans;
      } catch (const char *) {}
    }
  } else {
    for (float f = clb[nvar-1]; f < cub[nvar-1]; f = nextafterf(f, 1.0/0.0)) {
      vector<pair<expression *, pair<float, float> > > new_ineqs;
      FOR(i, ineqs.size()) new_ineqs.push_back(make_pair(
          constant_fold(subs(ineqs[i].first, nvar-1, f)),
          ineqs[i].second));
      try {
        vector<float> ans = dive(nvar-1, clb, cub, new_ineqs, preferred, 100);
        ans.push_back(f);
        return ans;
      } catch (const char *) {}
    }
  }
  throw ":(";
}

// Given a straight-line program `e`, an interval `[xlb, xub]` of `float`s, and
// a list `c` of coefficients, find at least one point `x` at which `e` with
// coefficients `c` yields an unacceptable function value.
int find_cuts(expression *e, float xlb, float xub, const vector<float> &c, vector<float> &testpoints) {
  float foo[c.size()+1];
  FOR(i, c.size()) foo[i+1] = c[i];
  int found = 0;
  FOR(i, 1000000) {
    float x = xlb + drand48() * (xub - xlb);
    foo[0] = x;
    float fx = eval(e, foo+1);
    float lb, ub;
    get_bounds(x, lb, ub);
    if (fx < lb || fx >= ub) testpoints.push_back(x), found++;
    if (found) return found;
  }
  for (float x = xub; x >= xlb && !found; x = nextafterf(x, -1.0/0.0)) {
    foo[0] = x;
    float fx = eval(e, foo+1);
    float lb, ub;
    get_bounds(x, lb, ub);
    if (fx < lb || fx >= ub) testpoints.push_back(x), found++;
  }
  return found;
}

// Driver for the heuristic.  Given a straight-line program `e`, the number of
// unspecified coefficients `nvar`, an interval `[xlb, xub]` of `float`s,
// bounds `[clb_0, cub_0] x ... x [clb_{nvar-1}, cub_{nvar-1}]` on the
// unspecified coefficients, and a nonempty vector `testpoints` of test points,
// try to find coefficients that yield an acceptable approximation.
int findit(expression *e, int nvar, float xlb, float xub,
    float *clb, float *cub, vector<float> &testpoints) {
  vector<float> coeffs(nvar);
  while (1) {
    try {
      lp_t lp(nvar, clb, cub);
      FOR(i, testpoints.size()) {
        gen_inequalities(e, nvar, clb, cub, testpoints[i], &lp);
      }
      FOR(i, nvar) {
        get_var_bounds(&lp, i, clb[i], cub[i]);
        cub[i] = nextafterf(cub[i], -1.0/0.0);
      }
      FOR(i, nvar) printf("%.6a <= c%i <= %.6a\n", clb[i], i, cub[i]);
    } catch (const char *) { return -2; } // infeasible.

    lp_t lp(nvar, clb, cub);
    vector<pair<expression *, pair<float, float> > > exprs;
    FOR(i, testpoints.size()) {
      float lower, upper;
      get_bounds(testpoints[i], lower, upper);
      exprs.push_back(make_pair(constant_fold(subs(e, -1, testpoints[i])),
                                make_pair(lower, upper)));
    }
    try {
      coeffs = dive(nvar, clb, cub, exprs, coeffs, 50);
    } catch (const char *) { return -1; } // infeasible or too hard to solve.
    printf("usces: %lli\n", get_cpu_usecs() - program_start);
    FOR(i, coeffs.size()) printf("float c%i = %a\n", i, coeffs[i]);
    if (!find_cuts(e, xlb, xub, coeffs, testpoints)) return 0;
    printf("vector<float> testpoints = {");
    FOR(i, testpoints.size()) printf(i?",%a":"%a", testpoints[i]);
    printf("};\n");
  }
}

int main() {
  program_start = get_cpu_usecs();
  auto fma = expression::fma;
  expression *x = expression::var(-1);
  expression *c3 = expression::var(6);
  expression *c5 = expression::var(5);
  expression *c7 = expression::var(4);
  expression *c9 = expression::var(3);
  expression *c11 = expression::var(2);
  expression *c13 = expression::var(1);
  expression *c15 = expression::var(0);

  expression *s = fma(x, x, expression::con(0));
  expression *tan_poly = fma(s, c15, c13);
  tan_poly = fma(s, tan_poly, c11);
  tan_poly = fma(s, tan_poly, c9);
  tan_poly = fma(s, tan_poly, c7);
  tan_poly = fma(s, tan_poly, c5);
  tan_poly = fma(s, tan_poly, c3);
  tan_poly = fma(s, tan_poly, expression::con(0));
  tan_poly = fma(x, tan_poly, x);

  float clb[8] = {-1,-1,-1,-1,-1,-1,-1,-1};
  float cub[8] = {1,1,1,1,1,1,1,1};

  vector<float> testpoints = {0x1p-1};

  int k = findit(tan_poly, 7, 1e-4, M_PI/4, clb, cub, testpoints);

  printf("%i\n", k);
}
