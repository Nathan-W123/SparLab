#include "sparlab/topopt/Mma.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <sstream>

namespace sparlab {

void MmaOptions::validate() const {
  if (!(move_limit > 0.0 && move_limit <= 1.0)) {
    throw ConfigError("mma.move_limit must lie in (0, 1]");
  }
  if (!(asymptote_init > 0.0 && asymptote_init <= 10.0)) {
    throw ConfigError("mma.asymptote_init must lie in (0, 10]");
  }
  if (!(asymptote_increase > 1.0) || !(asymptote_decrease > 0.0 && asymptote_decrease < 1.0)) {
    throw ConfigError(
        "mma.asymptote_increase must exceed 1 and mma.asymptote_decrease must lie in (0, 1)");
  }
  if (!(albefa > 0.0 && albefa < 1.0)) throw ConfigError("mma.albefa must lie in (0, 1)");
  if (!(raa0 >= 0.0)) throw ConfigError("mma.raa0 must be non-negative");
  if (!(c > 0.0) || !(d >= 0.0) || !(a0 > 0.0)) {
    throw ConfigError("mma.c and mma.a0 must be positive and mma.d non-negative");
  }
  if (!(epsimin > 0.0 && epsimin < 1.0)) throw ConfigError("mma.epsimin must lie in (0, 1)");
  if (!(constraint_scale_cap >= 1.0)) {
    throw ConfigError("mma.constraint_scale_cap must be at least 1");
  }
  if (max_inner_iterations < 10) throw ConfigError("mma.max_inner_iterations must be >= 10");
}

MmaOptimizer::MmaOptimizer(Index num_variables, Index num_constraints, const Vector& lower,
                           const Vector& upper, MmaOptions options)
    : n_(num_variables), m_(num_constraints), xmin_(lower), xmax_(upper), options_(options) {
  options_.validate();
  if (n_ < 1 || m_ < 1) throw ConfigError("MMA needs at least one variable and one constraint");
  if (xmin_.size() != n_ || xmax_.size() != n_) {
    throw ConfigError("MMA bound vectors must have one entry per variable");
  }
  for (Index j = 0; j < n_; ++j) {
    if (!(xmax_(j) > xmin_(j))) {
      std::ostringstream os;
      os << "MMA variable " << j << " has an empty box [" << xmin_(j) << ", " << xmax_(j)
         << "]; fixed variables must be removed before calling MMA";
      throw ConfigError(os.str());
    }
  }
  xold1_ = Vector::Zero(n_);
  xold2_ = Vector::Zero(n_);
  low_ = xmin_;
  upp_ = xmax_;
}

MmaStep MmaOptimizer::update(const Vector& x, Scalar f0, const Vector& df0dx,
                             const Vector& fval, const Matrix& dfdx) {
  (void)f0;  // only the gradient enters the separable approximation
  if (x.size() != n_ || df0dx.size() != n_ || fval.size() != m_ || dfdx.rows() != m_ ||
      dfdx.cols() != n_) {
    throw ConfigError("MMA received vectors of inconsistent size");
  }
  if (!x.allFinite() || !df0dx.allFinite() || !fval.allFinite() || !dfdx.allFinite()) {
    throw SolverError("MMA received non-finite function values or gradients");
  }
  ++iteration_;

  const Vector xmami_raw = xmax_ - xmin_;
  // Asymptotes: fixed distance for the first two iterations, then adapted to
  // the sign pattern of the last two moves.
  if (iteration_ <= 2) {
    low_ = x - options_.asymptote_init * xmami_raw;
    upp_ = x + options_.asymptote_init * xmami_raw;
  } else {
    for (Index j = 0; j < n_; ++j) {
      const Scalar zzz = (x(j) - xold1_(j)) * (xold1_(j) - xold2_(j));
      Scalar factor = 1.0;
      if (zzz > 0.0) factor = options_.asymptote_increase;
      if (zzz < 0.0) factor = options_.asymptote_decrease;
      low_(j) = x(j) - factor * (xold1_(j) - low_(j));
      upp_(j) = x(j) + factor * (upp_(j) - xold1_(j));
      const Scalar lowmin = x(j) - 10.0 * xmami_raw(j);
      const Scalar lowmax = x(j) - 0.01 * xmami_raw(j);
      const Scalar uppmin = x(j) + 0.01 * xmami_raw(j);
      const Scalar uppmax = x(j) + 10.0 * xmami_raw(j);
      low_(j) = std::clamp(low_(j), lowmin, lowmax);
      upp_(j) = std::clamp(upp_(j), uppmin, uppmax);
    }
  }

  // Bounds of the subproblem: inside the asymptotes with a margin, inside
  // the move limits, inside the box.
  Vector alfa(n_);
  Vector beta(n_);
  for (Index j = 0; j < n_; ++j) {
    alfa(j) = std::max({low_(j) + options_.albefa * (x(j) - low_(j)),
                        x(j) - options_.move_limit * xmami_raw(j), xmin_(j)});
    beta(j) = std::min({upp_(j) - options_.albefa * (upp_(j) - x(j)),
                        x(j) + options_.move_limit * xmami_raw(j), xmax_(j)});
  }

  // Approximation coefficients.
  const Vector xmami = xmami_raw.cwiseMax(1.0e-5);
  const Vector xmamiinv = xmami.cwiseInverse();
  const Vector ux1 = upp_ - x;
  const Vector xl1 = x - low_;
  const Vector ux2 = ux1.cwiseProduct(ux1);
  const Vector xl2 = xl1.cwiseProduct(xl1);
  const Vector uxinv = ux1.cwiseInverse();
  const Vector xlinv = xl1.cwiseInverse();

  Vector p0 = df0dx.cwiseMax(0.0);
  Vector q0 = (-df0dx).cwiseMax(0.0);
  const Vector pq0 = 0.001 * (p0 + q0) + options_.raa0 * xmamiinv;
  p0 = (p0 + pq0).cwiseProduct(ux2);
  q0 = (q0 + pq0).cwiseProduct(xl2);

  // Per-constraint scaling of badly violated rows (see MmaOptions).
  Vector row_scale = Vector::Ones(m_);
  const Scalar move = options_.move_limit * xmami_raw.maxCoeff();
  for (Index i = 0; i < m_; ++i) {
    const Scalar magnitude =
        std::max(std::abs(fval(i)), dfdx.row(i).cwiseAbs().maxCoeff() * move);
    if (magnitude > options_.constraint_scale_cap) {
      row_scale(i) = options_.constraint_scale_cap / magnitude;
    }
  }
  const Vector fscaled = row_scale.cwiseProduct(fval);
  const Matrix dscaled = row_scale.asDiagonal() * dfdx;

  Matrix p = dscaled.cwiseMax(0.0);
  Matrix q = (-dscaled).cwiseMax(0.0);
  const Matrix pq = 0.001 * (p + q) + options_.raa0 * Vector::Ones(m_) * xmamiinv.transpose();
  p = (p + pq) * ux2.asDiagonal();
  q = (q + pq) * xl2.asDiagonal();
  const Vector b = p * uxinv + q * xlinv - fscaled;

  MmaStep step = solve_subproblem(alfa, beta, p0, q0, p, q, b);
  step.max_change = (step.x - x).cwiseAbs().maxCoeff();
  // Multipliers and slacks back in the caller's scaling: s f_i <= 0 with
  // multiplier l_s is f_i <= 0 with multiplier s l_s.
  step.lambda = step.lambda.cwiseProduct(row_scale);
  step.y = step.y.cwiseQuotient(row_scale);
  step.constraint_scale = row_scale;

  xold2_ = xold1_;
  xold1_ = x;
  return step;
}

MmaStep MmaOptimizer::solve_subproblem(const Vector& alfa, const Vector& beta,
                                       const Vector& p0, const Vector& q0, const Matrix& p,
                                       const Matrix& q, const Vector& b) const {
  const Index n = n_;
  const Index m = m_;
  const Vector a = Vector::Zero(m);
  const Vector c = Vector::Constant(m, options_.c);
  const Vector d = Vector::Constant(m, options_.d);
  const Scalar a0 = options_.a0;

  // Primal-dual variables.
  Scalar epsi = 1.0;
  Vector x = 0.5 * (alfa + beta);
  Vector y = Vector::Ones(m);
  Scalar z = 1.0;
  Vector lam = Vector::Ones(m);
  Vector xsi = (x - alfa).cwiseInverse().cwiseMax(1.0);
  Vector eta = (beta - x).cwiseInverse().cwiseMax(1.0);
  Vector mu = (0.5 * c).cwiseMax(1.0);
  Scalar zet = 1.0;
  Vector s = Vector::Ones(m);

  // Residual of the relaxed KKT system for the current iterate.
  struct Residual {
    Scalar norm = 0.0;
    Scalar max = 0.0;
  };
  const auto residual = [&](const Vector& xv, const Vector& yv, Scalar zv, const Vector& lamv,
                            const Vector& xsiv, const Vector& etav, const Vector& muv,
                            Scalar zetv, const Vector& sv) {
    const Vector ux1 = upp_ - xv;
    const Vector xl1 = xv - low_;
    const Vector ux2 = ux1.cwiseProduct(ux1);
    const Vector xl2 = xl1.cwiseProduct(xl1);
    const Vector plam = p0 + p.transpose() * lamv;
    const Vector qlam = q0 + q.transpose() * lamv;
    const Vector gvec = p * ux1.cwiseInverse() + q * xl1.cwiseInverse();
    const Vector dpsidx = plam.cwiseQuotient(ux2) - qlam.cwiseQuotient(xl2);

    const Vector rex = dpsidx - xsiv + etav;
    const Vector rey = c + d.cwiseProduct(yv) - muv - lamv;
    const Scalar rez = a0 - zetv - a.dot(lamv);
    const Vector relam = gvec - a * zv - yv + sv - b;
    const Vector rexsi = xsiv.cwiseProduct(xv - alfa).array() - epsi;
    const Vector reeta = etav.cwiseProduct(beta - xv).array() - epsi;
    const Vector remu = muv.cwiseProduct(yv).array() - epsi;
    const Scalar rezet = zetv * zv - epsi;
    const Vector res = lamv.cwiseProduct(sv).array() - epsi;

    Residual r;
    Scalar sq = rex.squaredNorm() + rey.squaredNorm() + rez * rez + relam.squaredNorm() +
                rexsi.squaredNorm() + reeta.squaredNorm() + remu.squaredNorm() +
                rezet * rezet + res.squaredNorm();
    r.norm = std::sqrt(sq);
    r.max = std::max({rex.cwiseAbs().maxCoeff(), rey.cwiseAbs().maxCoeff(), std::abs(rez),
                      relam.cwiseAbs().maxCoeff(), rexsi.cwiseAbs().maxCoeff(),
                      reeta.cwiseAbs().maxCoeff(), remu.cwiseAbs().maxCoeff(),
                      std::abs(rezet), res.cwiseAbs().maxCoeff()});
    return r;
  };

  int total_iterations = 0;
  bool converged = true;
  Residual current;

  while (epsi > options_.epsimin) {
    current = residual(x, y, z, lam, xsi, eta, mu, zet, s);
    int ittt = 0;
    while (current.max > 0.9 * epsi && ittt < options_.max_inner_iterations) {
      ++ittt;
      ++total_iterations;

      const Vector ux1 = upp_ - x;
      const Vector xl1 = x - low_;
      const Vector ux2 = ux1.cwiseProduct(ux1);
      const Vector xl2 = xl1.cwiseProduct(xl1);
      const Vector ux3 = ux1.cwiseProduct(ux2);
      const Vector xl3 = xl1.cwiseProduct(xl2);
      const Vector uxinv1 = ux1.cwiseInverse();
      const Vector xlinv1 = xl1.cwiseInverse();
      const Vector uxinv2 = ux2.cwiseInverse();
      const Vector xlinv2 = xl2.cwiseInverse();
      const Vector plam = p0 + p.transpose() * lam;
      const Vector qlam = q0 + q.transpose() * lam;
      const Vector gvec = p * uxinv1 + q * xlinv1;
      const Matrix gg = p * uxinv2.asDiagonal() - q * xlinv2.asDiagonal();  // m x n
      const Vector dpsidx = plam.cwiseQuotient(ux2) - qlam.cwiseQuotient(xl2);

      const Vector xa = x - alfa;
      const Vector bx = beta - x;
      const Vector delx = dpsidx - epsi * xa.cwiseInverse() + epsi * bx.cwiseInverse();
      const Vector dely = c + d.cwiseProduct(y) - lam - epsi * y.cwiseInverse();
      const Scalar delz = a0 - a.dot(lam) - epsi / z;
      const Vector dellam = gvec - a * z - y - b + epsi * lam.cwiseInverse();

      Vector diagx = plam.cwiseQuotient(ux3) + qlam.cwiseQuotient(xl3);
      diagx = 2.0 * diagx + xsi.cwiseQuotient(xa) + eta.cwiseQuotient(bx);
      const Vector diagxinv = diagx.cwiseInverse();
      const Vector diagy = d + mu.cwiseQuotient(y);
      const Vector diagyinv = diagy.cwiseInverse();
      const Vector diaglam = s.cwiseQuotient(lam);
      const Vector diaglamyi = diaglam + diagyinv;

      Vector dx(n);
      Vector dlam(m);
      Scalar dz = 0.0;
      if (m < n) {
        // Reduce to the (lambda, z) system.
        const Vector blam = dellam + dely.cwiseQuotient(diagy) - gg * delx.cwiseQuotient(diagx);
        Vector bb(m + 1);
        bb.head(m) = blam;
        bb(m) = delz;
        Matrix alam = Matrix(diaglamyi.asDiagonal()) + gg * diagxinv.asDiagonal() * gg.transpose();
        Matrix aa = Matrix::Zero(m + 1, m + 1);
        aa.topLeftCorner(m, m) = alam;
        aa.topRightCorner(m, 1) = a;
        aa.bottomLeftCorner(1, m) = a.transpose();
        aa(m, m) = -zet / z;
        const Vector solut = aa.fullPivLu().solve(bb);
        dlam = solut.head(m);
        dz = solut(m);
        dx = -delx.cwiseQuotient(diagx) - (gg.transpose() * dlam).cwiseQuotient(diagx);
      } else {
        // Reduce to the (x, z) system (only for tiny problems).
        const Vector diaglamyiinv = diaglamyi.cwiseInverse();
        const Vector dellamyi = dellam + dely.cwiseQuotient(diagy);
        Matrix axx = Matrix(diagx.asDiagonal()) + gg.transpose() * diaglamyiinv.asDiagonal() * gg;
        const Scalar azz = zet / z + a.dot(a.cwiseQuotient(diaglamyi));
        const Vector axz = -gg.transpose() * a.cwiseQuotient(diaglamyi);
        const Vector bxv = delx + gg.transpose() * dellamyi.cwiseQuotient(diaglamyi);
        const Scalar bz = delz - a.dot(dellamyi.cwiseQuotient(diaglamyi));
        Matrix aa = Matrix::Zero(n + 1, n + 1);
        aa.topLeftCorner(n, n) = axx;
        aa.topRightCorner(n, 1) = axz;
        aa.bottomLeftCorner(1, n) = axz.transpose();
        aa(n, n) = azz;
        Vector bb(n + 1);
        bb.head(n) = -bxv;
        bb(n) = -bz;
        const Vector solut = aa.fullPivLu().solve(bb);
        dx = solut.head(n);
        dz = solut(n);
        dlam = (gg * dx).cwiseQuotient(diaglamyi) - dz * a.cwiseQuotient(diaglamyi) +
               dellamyi.cwiseQuotient(diaglamyi);
      }

      const Vector dy = -dely.cwiseQuotient(diagy) + dlam.cwiseQuotient(diagy);
      const Vector dxsi = -xsi + epsi * xa.cwiseInverse() - xsi.cwiseProduct(dx).cwiseQuotient(xa);
      const Vector deta = -eta + epsi * bx.cwiseInverse() + eta.cwiseProduct(dx).cwiseQuotient(bx);
      const Vector dmu = -mu + epsi * y.cwiseInverse() - mu.cwiseProduct(dy).cwiseQuotient(y);
      const Scalar dzet = -zet + epsi / z - zet * dz / z;
      const Vector ds = -s + epsi * lam.cwiseInverse() - s.cwiseProduct(dlam).cwiseQuotient(lam);

      // Step length keeping every positive variable positive.
      Scalar stmxx = 0.0;
      const auto ratio_max = [&](const Vector& dv, const Vector& v) {
        Scalar worst = -std::numeric_limits<Scalar>::max();
        for (Eigen::Index i = 0; i < v.size(); ++i) worst = std::max(worst, -1.01 * dv(i) / v(i));
        return worst;
      };
      stmxx = std::max({ratio_max(dy, y), -1.01 * dz / z, ratio_max(dlam, lam),
                        ratio_max(dxsi, xsi), ratio_max(deta, eta), ratio_max(dmu, mu),
                        -1.01 * dzet / zet, ratio_max(ds, s)});
      const Scalar stmalfa = ratio_max(dx, xa);
      Scalar stmbeta = -std::numeric_limits<Scalar>::max();
      for (Index j = 0; j < n; ++j) stmbeta = std::max(stmbeta, 1.01 * dx(j) / bx(j));
      const Scalar stminv = std::max({stmxx, stmalfa, stmbeta, 1.0});
      Scalar steg = 1.0 / stminv;

      const Vector xold = x;
      const Vector yold = y;
      const Scalar zold = z;
      const Vector lamold = lam;
      const Vector xsiold = xsi;
      const Vector etaold = eta;
      const Vector muold = mu;
      const Scalar zetold = zet;
      const Vector sold = s;

      // Backtracking on the residual norm.
      int itto = 0;
      Residual trial;
      trial.norm = 2.0 * current.norm;
      while (trial.norm > current.norm && itto < 50) {
        ++itto;
        x = xold + steg * dx;
        y = yold + steg * dy;
        z = zold + steg * dz;
        lam = lamold + steg * dlam;
        xsi = xsiold + steg * dxsi;
        eta = etaold + steg * deta;
        mu = muold + steg * dmu;
        zet = zetold + steg * dzet;
        s = sold + steg * ds;
        trial = residual(x, y, z, lam, xsi, eta, mu, zet, s);
        steg *= 0.5;
      }
      current = trial;
    }
    if (ittt >= options_.max_inner_iterations) {
      converged = false;
      log::warn("MMA subproblem: ", ittt, " Newton iterations at barrier ", epsi,
                " without reaching the residual target (", current.max, " > ",
                0.9 * epsi, ")");
    }
    epsi *= 0.1;
  }

  MmaStep step;
  step.x = x;
  step.lambda = lam;
  step.y = y;
  step.z = z;
  step.subproblem_iterations = total_iterations;
  step.subproblem_residual = current.max;
  step.subproblem_converged = converged;
  if (!converged) {
    std::ostringstream os;
    os << "the MMA subproblem solver did not converge (final residual " << current.max
       << " at barrier " << options_.epsimin
       << "); the approximation is badly scaled - check the objective and constraint "
          "scaling or loosen mma.epsimin";
    throw ConvergenceError(os.str());
  }
  return step;
}

}  // namespace sparlab
