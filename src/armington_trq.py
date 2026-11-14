import numpy as np
import matplotlib.pyplot as plt


GREEN, RED, TOTAL_COLOR = "#3f7d58", "#b1443a", "#1f1f1f"
LABEL_OFFSET = 0.26  # gap between a bar's end and its value label
YLIM_PAD = 1.4       # extra vertical room so labels don't clip
PLOT_STYLE = {
    "font.family": "serif",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.spines.left": False,
    "axes.grid": True,
    "axes.grid.axis": "y",
    "grid.color": "#e6e6e6",
    "grid.linewidth": 0.9,
    "axes.edgecolor": "#888888",
    "text.color": "#1a1a1a",
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
}


class PathMCP:
    """PATH solver (Ferris & Munson, 2000): non-monotone stabilization scheme
    for Mixed Complementarity Problems (MCP).

    Reformulates the MCP via the natural map (equivalent to Fischer-Burmeister,
    correctly handles boxes [lower, upper], including the degenerate case
    lower = upper):
        Phi(z) = z - proj_{[lower, upper]}(z - F(z)).
    Solves with semismooth Newton, a non-monotone line search and a
    Levenberg-Marquardt fallback (stabilization).
    """

    MEMORY = 10  # window of past merit values for the non-monotone line search

    def __init__(self, F, lower, upper, merit_tol=1e-10, maxiter=200):
        self.F = F
        self.lower = np.asarray(lower, dtype=float)
        self.upper = np.asarray(upper, dtype=float)
        self.merit_tol = merit_tol
        self.maxiter = maxiter

    def phi(self, z, Fz=None):
        if Fz is None:
            Fz = self.F(z)
        return z - np.clip(z - Fz, self.lower, self.upper)

    def merit(self, z, Fz=None):
        Phi = self.phi(z, Fz)
        return 0.5 * np.dot(Phi, Phi)

    def jac_F(self, z, rel=1e-6, absstep=1e-8):
        """Central finite-difference Jacobian of F."""
        m = len(z)
        J = np.zeros((m, m))
        for j in range(m):
            h = max(absstep, rel * abs(z[j]))
            zp = z.copy()
            zp[j] += h
            zm = z.copy()
            zm[j] -= h
            J[:, j] = (self.F(zp) - self.F(zm)) / (2 * h)
        return J

    def jac_phi(self, z, Fz):
        """Jacobian of Phi: rows of J_F where the projection is inactive,
        identity rows where it is clipped at a bound."""
        y = z - Fz
        interior = (y > self.lower) & (y < self.upper)
        return np.where(interior[:, None], self.jac_F(z), np.eye(len(z)))

    @staticmethod
    def _safe_solve(A, b):
        try:
            return np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            return np.linalg.lstsq(A, b, rcond=None)[0]

    def _linesearch(self, z, d, gdot, merit_ref):
        """Backtracking (Armijo) search against a reference merit value."""
        alpha = 1.0
        while alpha > 1e-14:
            z_trial = z + alpha * d
            if self.merit(z_trial) <= merit_ref + 1e-4 * alpha * gdot:
                return True, z_trial
            alpha *= 0.5
        return False, z

    def _levenberg_marquardt(self, z, J, Phi, merit_ref, lam):
        """Damped Gauss-Newton step on the merit function.
        Returns (ok, z, lam)."""
        grad = J.T @ Phi
        for _ in range(40):
            d = self._safe_solve(J.T @ J + lam * np.eye(len(z)), -grad)
            ok, z_trial = self._linesearch(z, d, np.dot(grad, d), merit_ref)
            if ok:
                return True, z_trial, lam
            lam *= 10.0
        return False, z, lam

    def _result(self, z, it, status, merit=None):
        if merit is None:
            merit = self.merit(z)
        return z, {"iter": it, "merit": merit, "status": status}

    def solve(self, z0):
        z = np.asarray(z0, dtype=float)
        merit_hist = []
        lam = 1e-8

        for it in range(self.maxiter):
            Fz = self.F(z)
            Phi = self.phi(z, Fz)
            merit = 0.5 * np.dot(Phi, Phi)
            merit_hist.append(merit)
            if merit <= self.merit_tol:
                return self._result(z, it, "converged", merit)

            merit_ref = max(merit_hist[-self.MEMORY:])

            # 1) semismooth Newton step
            J = self.jac_phi(z, Fz)
            d = self._safe_solve(J, -Phi)
            ok, z = self._linesearch(z, d, np.dot(Phi, J @ d), merit_ref)
            if ok:
                lam = max(1e-8, lam * 0.5)
                continue

            # 2) Levenberg-Marquardt stabilization
            ok, z, lam = self._levenberg_marquardt(
                z, J, Phi, merit_ref, max(lam, 1e-6))
            if ok:
                continue

            # 3) small steepest-descent step, or give up if it doesn't help
            z_trial = z - 1e-4 * (J.T @ Phi)
            if self.merit(z_trial) >= merit:
                return self._result(z, it, "stalled", merit)
            z = z_trial

        return self._result(z, self.maxiter, "maxiter")


class ArmingtonTRQ:
    """Armington model with a tariff-rate quota (TRQ), solved as an MCP.

    Total import demand is CES (Armington), with elasticity of substitution
    `sigma` across origins and aggregate demand elasticity `eta`. Each
    origin's supply responds to the received price with elasticity `eps[i]`.
    Regions with `qbar[i] > 0` face a quantitative quota: they pay `t_in[i]`
    in-quota and `t_out[i]` out-of-quota; the quota shadow price (omega) is
    endogenous and is determined in `solve()`.

    Parameters (vectors with one entry per region, except eta/sigma):
        regions : region names.
        v0      : base values; base quantities q0 = v0 (p0 = 1).
        eta     : total demand elasticity (scalar).
        sigma   : Armington elasticity of substitution (scalar).
        eps     : supply price elasticity by region.
        t0      : initial reference tariffs (ad valorem).
        qbar    : quota volume (0 = region without quota).
        t_in    : in-quota tariff (ad valorem).
        t_out   : out-of-quota tariff (ad valorem).

    Example:
        >>> model = ArmingtonTRQ(qbar=(0, 30, 0), t_out=(0, 0.8, 0))
        >>> model.solve().report()
    """

    def __init__(
        self,
        regions=("dom", "sub", "nsub"),
        v0=(60.0, 30.0, 10.0),
        eta=1.0,
        sigma=4.0,
        eps=(1.0, 10.0, 10.0),
        t0=(0.0, 0.0, 0.0),
        qbar=(0.0, 15.0, 0.0),
        t_in=(0.0, 0.1, 0.0),
        t_out=(0.0, 0.5, 0.0),
    ):
        self.regions = list(regions)
        n = len(self.regions)

        def as_vector(name, val):
            arr = np.asarray(val, dtype=float)
            if arr.shape != (n,):
                raise ValueError(
                    f"'{name}' must have {n} entries (one per region); "
                    f"got shape {arr.shape}."
                )
            return arr

        self.v0 = as_vector("v0", v0)
        self.eps = as_vector("eps", eps)
        self.t0 = as_vector("t0", t0)
        self.qbar = as_vector("qbar", qbar)
        self.t_in = as_vector("t_in", t_in)
        self.t_out = as_vector("t_out", t_out)
        self.eta = float(eta)
        self.sigma = float(sigma)

        self.q0 = self.v0.copy()
        self.share = self.v0 / self.v0.sum()
        self.quota_idx = np.where(self.qbar > 0)[0]

        if np.any(self.t_out[self.quota_idx] < self.t_in[self.quota_idx]):
            raise ValueError("'t_out' must be >= 't_in' for regions with a quota.")

        # Filled in by solve()
        self.p_eq = None
        self.omega = None
        self.info = None
        self.vars_out = None
        self.qd_out = None
        self.price_change_pct = None
        self.quantity_change_pct = None

    def _check_solved(self):
        if self.p_eq is None:
            raise RuntimeError("Call solve() before this method.")

    def compute_variables(self, p, omega):
        """Prices, demand and supply given prices `p` and quota rents `omega`
        (full-length vectors; omega is zero for regions without a quota)."""
        p = np.asarray(p, dtype=float)
        t_eff = self.t_in + omega  # omega = 0 where there is no quota
        taup = (1 + t_eff) / (1 + self.t0)
        P = (self.share * p ** (1 - self.sigma)).sum() ** (1 / (1 - self.sigma))
        Q = P ** (-abs(self.eta))
        return {
            "P": P,
            "qd": (p / P) ** (-self.sigma) * Q,
            "qs": (p / taup) ** self.eps,
            "t_eff": t_eff,
        }

    def solve(self, merit_tol=1e-10, maxiter=200):
        """Solve the equilibrium; store results as attributes."""
        n = len(self.regions)
        quota_idx = self.quota_idx
        n_quota = len(quota_idx)

        def expand_omega(omega_quota):
            omega = np.zeros(n)
            omega[quota_idx] = omega_quota
            return omega

        # Unknowns z = (prices, quota rents); F = (excess demand, quota gap)
        def F(z):
            out = self.compute_variables(z[:n], expand_omega(z[n:]))
            excess_demand = out["qd"] - out["qs"]
            quota_gap = (self.qbar[quota_idx] / self.q0[quota_idx]
                         - out["qd"][quota_idx])
            return np.concatenate([excess_demand, quota_gap])

        lower = np.concatenate([np.full(n, -np.inf), np.zeros(n_quota)])
        upper = np.concatenate([np.full(n, np.inf),
                                self.t_out[quota_idx] - self.t_in[quota_idx]])
        z0 = np.concatenate([np.ones(n), np.zeros(n_quota)])
        solver = PathMCP(F, lower, upper, merit_tol=merit_tol, maxiter=maxiter)
        z, self.info = solver.solve(z0)

        self.p_eq = z[:n]
        self.omega = expand_omega(z[n:])
        self.vars_out = self.compute_variables(self.p_eq, self.omega)
        self.qd_out = self.vars_out["qd"] * self.q0
        self.price_change_pct = (self.p_eq - 1.0) * 100.0
        self.quantity_change_pct = (self.qd_out - self.q0) / self.q0 * 100.0
        return self

    def welfare(self):
        """Welfare decomposition for the domestic origin (region 0)."""
        self._check_solved()
        dom = 0
        eps_dom = self.eps[dom]
        base_exp = self.q0.sum()

        # Delta CS (total demand curve / expenditure function E = u * P)
        dCS = -base_exp * (self.vars_out["P"] - 1.0)

        # Delta PS (domestic production; taup = 1 since t0[0] = 0)
        dPS = self.q0[dom] * (self.p_eq[dom] ** (eps_dom + 1.0) - 1.0) / (eps_dom + 1.0)

        # Tariff revenue + quota rent (tariffed origins = regions with a quota)
        quota_imports = self.qd_out[self.quota_idx]
        tar_rev = np.sum(self.t_in[self.quota_idx] * quota_imports)
        quota_rent = np.sum(self.omega[self.quota_idx] * quota_imports)

        return {
            "dCS": dCS,
            "dPS": dPS,
            "tar_rev": tar_rev,
            "quota_rent": quota_rent,
            "total": dCS + dPS + tar_rev + quota_rent,
            "base_exp": base_exp,
        }

    def _regime(self, i):
        """Label describing how the quota behaves for region `i`."""
        if self.qbar[i] <= 0:
            return "no quota"
        tariff_gap = self.t_out[i] - self.t_in[i]
        if tariff_gap <= 1e-9:
            return "non-binding"
        if self.omega[i] <= 1e-6:
            return "in-quota"
        if self.omega[i] >= tariff_gap - 1e-6:
            return "out-of-quota"
        return "binding quota"

    def report(self):
        """Print the per-region table and the welfare decomposition."""
        self._check_solved()
        w = self.welfare()

        print("\nResults by region (PATH solver - non-monotone stabilization)")
        print("-" * 78)
        print(f"{'Region':<10}{'chg. p (%)':>12}{'chg. q (%)':>12}{'Quota':>8}"
              f"{'Eff. tariff':>12}{'Regime':>14}")
        print("-" * 78)
        for i, region in enumerate(self.regions):
            quota_str = f"{self.qbar[i]:g}" if self.qbar[i] > 0 else "-"
            print(f"{region:<10}{self.price_change_pct[i]:>11.2f}%"
                  f"{self.quantity_change_pct[i]:>11.2f}%{quota_str:>8}"
                  f"{self.vars_out['t_eff'][i]:>12.4f}{self._regime(i):>14}")
        print("-" * 78)
        print(f"PATH: status={self.info['status']}  iterations={self.info['iter']}  "
              f"merit={self.info['merit']:.3e}")

        def pct(x):
            return 100.0 * x / w["base_exp"]

        print("\n=== Welfare of the domestic origin (Brazil) - vs base p=1, t0=0 ===")
        print(f"ΔCS consumers         = {w['dCS']:8.4f}  ({pct(w['dCS']):6.2f}% base spend)")
        print(f"ΔPS producers (dom)   = {w['dPS']:8.4f}  ({pct(w['dPS']):6.2f}% base spend)")
        print(f"Tariff revenue        = {w['tar_rev']:8.4f}  (t_in x quota imports)")
        print(f"Quota rent            = {w['quota_rent']:8.4f}  (ω x quota imports)")
        print(f"Total welfare         = {w['total']:8.4f}  ({pct(w['total']):6.2f}% base spend)")

    def plot_welfare(self, save=None):
        """Waterfall of the welfare decomposition; saved to `save` if given."""
        w = self.welfare()
        labels = [
            "Consumers\n(ΔCS)",
            "Producers\n(ΔPS)",
            "Tariff\nrevenue",
            "Quota\nrent",
            "Total\nwelfare",
        ]
        # Everything is plotted as % of base spend
        to_pct = 100.0 / w["base_exp"]
        components = to_pct * np.array(
            [w["dCS"], w["dPS"], w["tar_rev"], w["quota_rent"]])
        total = to_pct * w["total"]
        tips = np.cumsum(components)          # running total after each bar
        bottoms = tips - components           # where each bar starts
        n_bars = len(labels)
        width = 0.62

        with plt.rc_context(PLOT_STYLE):
            fig, ax = plt.subplots(figsize=(9, 5.6))

            # Floating component bars, then the total bar from zero
            ax.bar(range(n_bars - 1), components, bottom=bottoms, width=width,
                   color=[RED if v < 0 else GREEN for v in components], zorder=3)
            ax.bar(n_bars - 1, total, width=width, color=TOTAL_COLOR,
                   zorder=3)

            # Value labels just beyond each bar's free end
            ends = np.append(tips, total)
            values = np.append(components, total)
            for i, (end, v) in enumerate(zip(ends, values)):
                ax.text(i, end + np.sign(v or 1) * LABEL_OFFSET,
                        f"{v:+.2f}%",
                        ha="center", va="bottom" if v >= 0 else "top",
                        fontsize=8, color="#1C1B1B")

            # Dashed connectors between consecutive bars
            for i, level in enumerate(tips):
                ax.plot([i + width / 2, i + 1 - width / 2], [level, level],
                        color="#999999", lw=1, ls="--", zorder=2)

            lo, hi = ax.get_ylim()
            ax.set_ylim(lo - YLIM_PAD, hi + YLIM_PAD)
            ax.axhline(0, color="#333333", lw=1, zorder=2)
            ax.set_xticks(range(n_bars))
            ax.set_xticklabels(labels, fontsize=9.5)
            ax.set_xlim(-0.6, n_bars - 0.4)
            ax.set_ylabel("Welfare change (% of base spend)", fontsize=10)
            ax.set_title("Welfare impact of the TRQ on the domestic origin",
                         fontsize=14, fontweight="bold", loc="left", pad=30,
                         color="#1a1a1a")
            ax.text(0, 1.005,
                    "Change in equivalent utility, as % of base spend (p = 1, t₀ = 0).",
                    transform=ax.transAxes, fontsize=9.5, color="#666666",
                    ha="left")
            fig.text(0.01, 0.02,
                     "Source: Armington-TRQ model (partial equilibrium, PATH solver).",
                     fontsize=8, color="#999999", style="italic")
            plt.tight_layout(rect=[0, 0.03, 1, 1])

        if save is not None:
            fig.savefig(save, dpi=150)
            plt.close(fig)
        else:
            plt.show()


if __name__ == "__main__":
    model = ArmingtonTRQ(
        t_in=[0.0, 0.1, 0.0],
        t_out=[0.0, 0.5, 0.0]
    )
    model.solve()
    model.report()
    model.plot_welfare()
