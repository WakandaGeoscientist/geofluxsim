# =============================================================================
#  GeoFluxSim with David — turbidity currents, their deposits, and the physics
#  Copy this whole file into ONE Jupyter cell and run it.
#  Requires: numpy, matplotlib, ipywidgets   (pip install numpy matplotlib ipywidgets)
#
#  Layout (top to bottom)
#    1. Run bar: Play/Pause, Reset, probe position, playback speed, live readouts
#    2. Parameters: sediment, fluid, flow and model coefficients (first thing you see)
#    3. Views: side view, 3D view, probe camera. Show the flow as Realistic,
#       Concentration or Velocity (colour-mapped, with the probe profile in white);
#       toggle the flow and the deposit layers; deposit coloured by composition or size
#    4. All plots in one dashboard: velocity & concentration profiles, front position,
#       sediment budget, deposit cross-section, core, downstream distribution,
#       live force/flux diagram, settling-velocity curves, and the equation table
#
#  Physics (layer-averaged, Parker et al. 1986), six grain classes
#  (siliciclastic & carbonate × fine / median / coarse from the sorting σg):
#    dh/dt + d(Uh)/dx = e_w|U|,   e_w = a / sqrt(1 + 718 Ri^2.4),  Ri = g'h/U²
#    d(Uh)/dt + d(U²h + ½g'h²)/dx = g'hS − Cd U|U|,   g' = g Σ R_i C_i
#    d(C_i h)/dt + d(U C_i h)/dx = − r0 w_s,i (1−C)^4.65 C_i
#  w_s,i: Ferguson & Church (2004) with grain density, shape (C1, C2), ambient
#  density and viscosity ν(T). Deposits stored in 1.5 s layers (stratigraphy).
# =============================================================================
import asyncio
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

G, KAPPA, HIND_N, ETA_M = 9.81, 0.41, 4.65, 0.2
H_MIN, EPS_H, CFL = 1.0e-4, 1.0e-3, 0.40
LAYER_DT, NR = 1.5, 128
BIN_W = np.array([0.25, 0.5, 0.25])
SH_SILI = {"Rounded quartz": (18.0, 1.0), "Angular quartz": (20.0, 1.2), "Platy, mica-rich": (24.0, 1.5)}
SH_CARB = {"Ooids / peloids (near-spherical)": (18.0, 0.6),
           "Bioclastic fragments (angular)": (20.0, 1.2),
           "Plates & rods (skeletal)": (24.0, 1.5)}
DEFAULTS = dict(L=8.0, W=0.6, nx=240, t_max=90.0,
                D_s=60.0, rho_s=2650.0, shape_s="Rounded quartz",
                D_c=60.0, rho_c=2300.0, shape_c="Bioclastic fragments (angular)",
                C0=5.0, fc=50.0, sort=1.6,
                T_w=20.0, rho_a=1000.0,
                slope=3.0, U0=0.20, h0=0.08, T_in=20.0,
                Cd=0.004, ea=0.075, r0=1.6)


def nu_water(T):
    """Kinematic viscosity of water (m²/s) at temperature T (°C)."""
    return 1.79e-6 / (1 + 0.03368 * T + 0.000221 * T * T)


def settling_velocity(D_um, rho_s, C1, C2, nu, rho_a):
    """Ferguson & Church (2004) settling velocity (m/s)."""
    D = np.asarray(D_um, float) * 1e-6
    R = (rho_s - rho_a) / rho_a
    if R <= 0:
        return np.zeros_like(D)
    return R * G * D**2 / (C1 * nu + np.sqrt(0.75 * C2 * R * G * D**3))


def stokes_velocity(D_um, rho_s, nu, rho_a):
    D = np.asarray(D_um, float) * 1e-6
    return max((rho_s - rho_a) / rho_a, 0) * G * D**2 / (18 * nu)


def profile_k(r0):
    """Exponent k of the concentration profile with near-bed ratio r0: k/(1−e^−k) = r0."""
    k = 1.0
    for _ in range(40):
        e = np.exp(-k)
        f = k / (1 - e) - r0
        d = (1 - e - k * e) / (1 - e) ** 2
        k = max(k - f / d, 1e-3)
    return k


def _u_raw(eta):
    eta = np.maximum(eta, 0.0)
    return np.where(eta <= ETA_M, (eta / ETA_M) ** (1 / 6), np.exp(-1.5 * ((eta - ETA_M) / (1 - ETA_M)) ** 2))

_e = (np.arange(4000) + 0.5) * 3 / 4000
U_NORM = float(np.sum(_u_raw(_e)) * 3 / 4000)

def u_shape(eta):
    return _u_raw(eta) / U_NORM            # ∫ u dz = U h


# ======================================================================= model
class TurbidityModel:
    def __init__(self, **params):
        self.p = dict(DEFAULTS)
        self.p.update(params)
        self.reset()

    def reset(self):
        p = self.p
        self.L, self.nx = float(p["L"]), int(p["nx"])
        self.dx = self.L / self.nx
        self.x = (np.arange(self.nx) + 0.5) * self.dx
        self.S = p["slope"] / 100.0
        self.zb = self.S * (self.L - self.x)
        self.nu = nu_water(p["T_w"])
        # six classes: siliciclastic fine/median/coarse, carbonate fine/median/coarse
        sg = p["sort"]
        self.comp = np.array([0, 0, 0, 1, 1, 1])
        self.D = np.array([p["D_s"] * sg ** (b - 1) for b in range(3)] + [p["D_c"] * sg ** (b - 1) for b in range(3)])
        self.rho = np.array([p["rho_s"]] * 3 + [p["rho_c"]] * 3)
        shapes = [SH_SILI[p["shape_s"]]] * 3 + [SH_CARB[p["shape_c"]]] * 3
        frac = np.concatenate([(1 - p["fc"] / 100) * BIN_W, p["fc"] / 100 * BIN_W])
        self.R = (self.rho - p["rho_a"]) / p["rho_a"]
        self.ws0 = np.array([float(settling_velocity(self.D[k], self.rho[k], *shapes[k], self.nu, p["rho_a"]))
                             for k in range(6)])
        self.Cin = p["C0"] / 100 * frac
        self.kp = profile_k(p["r0"])
        self.ap = self.kp / (1 - np.exp(-self.kp))
        self.h = np.full(self.nx, H_MIN)
        self.q = np.zeros(self.nx)
        self.m = np.zeros((6, self.nx))
        self.dep = np.zeros((6, self.nx))
        self.layers = []
        self.injected = self.exited = 0.0
        self.t = self.front = self.phase = 0.0
        self.t_hist, self.xf_hist = [0.0], [0.0]
        self.budget = [(0.0, 0.0, 0.0, 0.0)]
        self.done = False

    @staticmethod
    def vel(h, q):
        h2 = h * h
        return 2.0 * h * q / (h2 + np.maximum(h2, EPS_H**2))

    def step(self):
        p, h, q, m = self.p, self.h, self.q, self.m
        if self.t < p["T_in"]:
            hL, qL, mL = p["h0"], p["h0"] * p["U0"], p["h0"] * self.Cin
        else:
            hL, qL, mL = h[0], -q[0], m[:, 0]
        hE = np.concatenate(([hL], h, [h[-1]]))
        qE = np.concatenate(([qL], q, [max(q[-1], 0.0)]))
        mE = np.concatenate((np.asarray(mL)[:, None], m, m[:, -1:]), axis=1)
        UE = self.vel(hE, qE)
        bh = G * (self.R[:, None] * mE).sum(0)
        F1 = qE * UE + 0.5 * bh * hE
        c = np.abs(UE) + np.sqrt(np.maximum(bh, 0.0))
        a = np.maximum(c[:-1], c[1:])
        fh = 0.5 * (qE[:-1] + qE[1:]) - 0.5 * a * (hE[1:] - hE[:-1])
        fq = 0.5 * (F1[:-1] + F1[1:]) - 0.5 * a * (qE[1:] - qE[:-1])
        Fm = UE * mE
        fm = 0.5 * (Fm[:, :-1] + Fm[:, 1:]) - 0.5 * a * (mE[:, 1:] - mE[:, :-1])
        dt = min(CFL * self.dx / max(c.max(), 1e-3), 0.05)
        r = dt / self.dx
        h_n = np.maximum(h - r * (fh[1:] - fh[:-1]), H_MIN)
        q_n = q - r * (fq[1:] - fq[:-1])
        m_n = np.maximum(m - r * (fm[:, 1:] - fm[:, :-1]), 0.0)
        self.injected += fm[:, 0].sum() * dt
        self.exited += np.maximum(fm[:, -1], 0.0).sum() * dt
        U = self.vel(h_n, q_n)
        bh = G * (self.R[:, None] * m_n).sum(0)
        Ri = np.minimum(bh / np.maximum(U * U, 1e-10), 1e3)
        h_n = h_n + p["ea"] / np.sqrt(1 + 718 * Ri**2.4) * np.abs(U) * dt
        q_n = (q_n + bh * self.S * dt) / (1 + dt * p["Cd"] * np.abs(U) / h_n)
        Ct = np.clip(m_n.sum(0) / h_n, 0, 0.6)
        ws = self.ws0[:, None] * ((1 - Ct) ** HIND_N)[None, :]
        m_new = m_n * np.exp(-p["r0"] * ws * dt / h_n)
        d = m_n - m_new
        self.dep += d
        li = int(self.t // LAYER_DT)
        while len(self.layers) <= li:
            self.layers.append(np.zeros((6, self.nx)))
        self.layers[li] += d
        self.h, self.q, self.m = h_n, q_n, np.minimum(m_new, 0.5 * h_n)
        self.t += dt
        idx = np.nonzero(self.m.sum(0) > 1e-3 * self.Cin.sum() * p["h0"])[0]
        self.front = min(self.x[idx[-1]] + 0.5 * self.dx, self.L) if idx.size else 0.0
        if self.t >= p["t_max"]:
            self.done = True

    def advance(self, sim_seconds):
        t_end = self.t + sim_seconds
        while self.t < t_end and not self.done:
            self.step()
        self.t_hist.append(self.t)
        self.xf_hist.append(self.front)
        inside = self.m.sum(0) / self.h > 1e-3 * self.Cin.sum()
        self.phase += (self.vel(self.h, self.q)[inside].mean() if inside.any() else 0.15) * sim_seconds
        inj = max(self.injected, 1e-12)
        self.budget.append((self.t, 100 * self.m.sum() * self.dx / inj,
                            100 * self.dep.sum() * self.dx / inj, 100 * self.exited / inj))

    def probe(self, xp):
        h = float(np.interp(xp, self.x, self.h))
        U = float(self.vel(np.array(h), np.array(np.interp(xp, self.x, self.q))))
        C = np.array([np.interp(xp, self.x, self.m[k]) for k in range(6)]) / h
        return dict(h=h, U=U, C=C, Cs=C[:3].sum(), Cc=C[3:].sum(), Ct=C.sum())

    def front_speed(self):
        if len(self.t_hist) < 4:
            return 0.0
        t, x = np.array(self.t_hist[-4:]), np.array(self.xf_hist[-4:])
        return float(np.polyfit(t, x, 1)[0]) if np.ptp(t) > 0 else 0.0

    def physics(self, xp):
        """Value of every term of the governing equations at the probe."""
        p, pr = self.p, self.probe(xp)
        i = int(np.clip(round(xp / self.dx - 0.5), 1, self.nx - 2))
        B = float((self.R * pr["C"]).sum()); gp = G * B; h, U = pr["h"], pr["U"]
        P = 0.5 * G * (self.R[:, None] * self.m).sum(0) * self.h
        Ri = gp * h / max(U * U, 1e-10)
        ew = p["ea"] / np.sqrt(1 + 718 * min(Ri, 1e3) ** 2.4)
        hind = (1 - min(pr["Ct"], 0.6)) ** HIND_N
        dep_k = p["r0"] * self.ws0 * hind * pr["C"]
        us = np.sqrt(p["Cd"]) * abs(U)
        return dict(pr=pr, B=B, gp=gp, h=h, U=U, grav=gp * h * self.S, press=-(P[i + 1] - P[i - 1]) / (2 * self.dx),
                    fric=-p["Cd"] * U * abs(U), entr=-ew * abs(U) * U, ew=ew, Ri=Ri,
                    Fr=abs(U) / np.sqrt(max(gp * h, 1e-12)), Re=abs(U) * h / self.nu, hind=hind,
                    depS=dep_k[:3].sum(), depC=dep_k[3:].sum(),
                    rouse=[self.ws0[1] * hind / (KAPPA * max(us, 1e-6)), self.ws0[4] * hind / (KAPPA * max(us, 1e-6))],
                    active=h > 0.003 and pr["Ct"] > 1e-4)


# ================================================================ deposits
COL_S, COL_MID, COL_C = np.array([224, 165, 38]), np.array([238, 232, 216]), np.array([31, 111, 178])
FINE, COARSE = np.array([244, 234, 206]), np.array([112, 66, 28])


def deposit_colour(cf, lnD, rng, mode):
    """Colour layers; rng = (cf0, cf1, lnD0, lnD1) is the deposit's own range (colours are stretched to it)."""
    cf0, cf1, d0, d1 = rng
    if mode == "Composition":
        t = (np.asarray(cf) - cf0) / (cf1 - cf0) if cf1 > cf0 + 0.02 else np.asarray(cf)
        t = np.clip(t, 0, 1)[..., None]
        lo = COL_S + (COL_MID - COL_S) * np.clip(t * 2, 0, 1)
        hi = COL_MID + (COL_C - COL_MID) * np.clip(t * 2 - 1, 0, 1)
        return np.where(t < 0.5, lo, hi) / 255
    t = np.clip((np.asarray(lnD) - d0) / (d1 - d0), 0, 1)[..., None] if d1 > d0 + 1e-3 else 0.5
    return (FINE + (COARSE - FINE) * t) / 255


def build_section(mdl, porosity, mode):
    """Layered deposit: returns dict with thickness (mm), image (NR, nx, 4) bottom row first, etc."""
    mmf = 1e3 / (1 - porosity)
    lnD = np.log(mdl.D)
    lnDmin, lnDmax = lnD.min(), lnD.max()
    thick = mdl.dep.sum(0) * mmf
    out = dict(thick=thick, mmS=mdl.dep[:3].sum(0) * mmf, mmC=mdl.dep[3:].sum(0) * mmf,
               meanD=np.where(thick > 0, np.exp((mdl.dep * lnD[:, None]).sum(0) / np.maximum(mdl.dep.sum(0), 1e-30)), np.nan),
               dmax=max(thick.max(), 0.01), lnDmin=lnDmin, lnDmax=lnDmax)
    if not mdl.layers:
        out.update(img=np.zeros((NR, mdl.nx, 4)), top=np.tile([0.35, 0.36, 0.37], (mdl.nx, 1)), layers=None,
                   rng=(0, 1, lnDmin, lnDmax))
        return out
    Ls = np.array(mdl.layers)                                  # (NL, 6, nx)
    v = Ls.sum(1)                                              # (NL, nx)
    cf = np.where(v > 0, Ls[:, 3:].sum(1) / np.maximum(v, 1e-30), 0)
    lD = np.where(v > 0, (Ls * lnD[None, :, None]).sum(1) / np.maximum(v, 1e-30), lnD.mean())
    big = v > 1e-3 * v.max()
    rng = (cf[big].min(), cf[big].max(), lD[big].min(), lD[big].max()) if big.any() else (0, 1, lnDmin, lnDmax)
    out["rng"] = rng
    col = deposit_colour(cf, lD, rng, mode)                   # (NL, nx, 3)
    cum = np.cumsum(v * mmf, axis=0)                           # (NL, nx)
    hgt = (np.arange(NR) + 0.5) / NR * out["dmax"]
    idx = np.minimum((cum[:, None, :] < hgt[None, :, None]).sum(0), len(Ls) - 1)   # (NR, nx)
    img = col[idx, np.arange(mdl.nx)[None, :]]                 # (NR, nx, 3)
    row_h = out["dmax"] / NR
    lay_top = np.take_along_axis(cum, idx, axis=0)
    timeline = ((idx % 5) == 4) & (np.abs(hgt[:, None] - lay_top) < row_h * 0.9) & (lay_top < thick[None, :] - row_h)
    img = np.where(timeline[..., None], img * 0.62, img)
    inside = hgt[:, None] <= thick[None, :]
    last = np.where(v > 0, np.arange(len(Ls))[:, None], -1).max(0)
    top = np.where(last[:, None] >= 0, col[np.maximum(last, 0), np.arange(mdl.nx)], np.array([0.35, 0.36, 0.37]))
    out.update(img=np.dstack([img, inside.astype(float)]), top=top, layers=(v * mmf, cf, np.exp(lD)))
    return out


# ================================================================ renderer
def _periodic_noise(shape, beta=3.0, seed=0, kmax=None):
    rng = np.random.default_rng(seed)
    ky, kx = np.fft.fftfreq(shape[0])[:, None], np.fft.fftfreq(shape[1])[None, :]
    k = np.sqrt(kx**2 + ky**2); k[0, 0] = 1.0
    amp = k ** (-beta / 2); amp[0, 0] = 0.0
    if kmax is not None:
        amp = amp * np.exp(-(k / kmax) ** 2)
    n = np.real(np.fft.ifft2(amp * np.exp(2j * np.pi * rng.random(shape))))
    return np.clip(0.5 + 0.17 * (n - n.mean()) / n.std(), 0, 1)

NOISE = _periodic_noise((256, 1024), seed=7)
NOISE_LO = _periodic_noise((256, 1024), seed=7, kmax=0.05)
NOISE2 = _periodic_noise((64, 512), beta=2.6, seed=3)


def _bilinear(tex, u, v):
    H, W = tex.shape
    u0, v0 = np.floor(u).astype(int), np.floor(v).astype(int)
    fu, fv = u - u0, v - v0
    u0 %= W; v0 %= H
    u1, v1 = (u0 + 1) % W, (v0 + 1) % H
    return (tex[v0, u0] * (1 - fu) + tex[v0, u1] * fu) * (1 - fv) + (tex[v1, u0] * (1 - fu) + tex[v1, u1] * fu) * fv


WATER_TOP, WATER_DEEP = np.array([0.70, 0.87, 0.96]), np.array([0.26, 0.50, 0.66])
ALB_SILI, ALB_CARB = np.array([0.70, 0.63, 0.52]), np.array([0.80, 0.79, 0.74])
LIGHT, AMBIENT = np.array([1.30, 1.22, 1.08]), np.array([0.13, 0.14, 0.15])
SIG_VIEW, SIG_LIGHT = 150.0, 1200.0


def field_cmap(mode):
    return plt.get_cmap("inferno" if mode == "Concentration" else "viridis")


def field_value(mode, c, u, vmax_c, vmax_u):
    if mode == "Concentration":
        return np.clip(c / vmax_c, 0, 1)
    return np.clip(u / vmax_u, 0, 1)


def render_flow(mdl, sec, exag, x0, x1, z0, z1, nxp, nzp, return_alpha=False,
                mode="Realistic", show_flow=True, show_dep=True, vmax_c=0.08, vmax_u=0.4):
    """Side view of current + layered deposit (top row = z1).
    mode 'Realistic' draws the flow like the lab footage; 'Concentration' / 'Velocity'
    colour the flow by that field (as a slice through the flume's centre line)."""
    xs, zs = np.linspace(x0, x1, nxp), np.linspace(z1, z0, nzp)
    dz = (z1 - z0) / (nzp - 1)
    h = np.interp(xs, mdl.x, mdl.h)
    mk = np.array([np.interp(xs, mdl.x, mdl.m[k]) for k in range(6)])
    C = mk.sum(0) / h
    fc = np.where(mk.sum(0) > 0, mk[3:].sum(0) / np.maximum(mk.sum(0), 1e-12), 0.5)
    U = np.interp(xs, mdl.x, mdl.vel(mdl.h, mdl.q))
    zb = mdl.S * (mdl.L - xs)
    dmm = np.interp(xs, mdl.x, sec["thick"]) if show_dep else np.zeros_like(xs)
    dh = dmm * 1e-3 * exag
    zeta_b = zs[:, None] - zb[None, :]
    zeta = zeta_b - dh[None, :]
    s = zeta / h[None, :]
    active = (h > 0.004) & (C > 2e-4) & show_flow
    ok = (zeta >= 0) & active[None, :]
    smooth = np.where(ok, C * mdl.ap * np.exp(-mdl.kp * np.maximum(s, 0)) * 0.5 * (1 - np.tanh((s - 1) / 0.12)), 0.0)
    k = np.clip((zs - z0) / max(z1 - z0, 1e-9), 0, 1)[:, None, None]
    bg = np.broadcast_to(WATER_DEEP * (1 - k) + WATER_TOP * k, (nzp, nxp, 3)).copy()
    if mode == "Realistic":
        u = (xs - mdl.phase) / 1.6 * 1024
        v = (s - mdl.t * 0.06) / 5.0 * 256
        tex = NOISE_LO if (x1 - x0) / nxp / 1.6 * 1024 > 2.0 else NOISE
        u = u[None, :] + 40.0 * (_bilinear(NOISE_LO, u[None, :] * 0.35, v * 0.5 + 77) - 0.5)
        n = _bilinear(tex, u, v)
        ramp = np.clip((s - 0.05) / 0.85, 0, 1); ramp = ramp * ramp * (3 - 2 * ramp)
        se = np.maximum(s - 1.4 * (n - 0.5) * ramp, 0.0)
        dens = np.where(ok, C * mdl.ap * np.exp(-mdl.kp * se) * 0.5 * (1 - np.tanh((se - 1) / 0.07)) * (0.5 + n), 0.0)
        T_light = np.exp(-SIG_LIGHT * np.cumsum(smooth, axis=0) * dz)
        alpha = 1.0 - np.exp(-SIG_VIEW * dens * mdl.p["W"])
        alb = ALB_SILI * (1 - fc)[None, :, None] + ALB_CARB * fc[None, :, None]
        lit = alb * (LIGHT * T_light[..., None] + AMBIENT)
    else:
        uf = np.where(ok, U[None, :] * u_shape(np.maximum(s, 0)), 0.0)
        val = field_value(mode, smooth, uf, vmax_c, vmax_u)
        lit = field_cmap(mode)(0.12 + 0.86 * val if mode == "Concentration" else val)[..., :3]
        g = np.abs(np.mod(val * 8, 1) - 0.5)
        lit = lit * np.where(g > 0.44, 0.7, 1.0)[..., None]
        alpha = np.clip(smooth / (0.02 * vmax_c), 0, 1)
    in_dep = (zeta_b >= 0) & (zeta_b < dh[None, :])
    row = np.clip((zeta_b / max(exag, 1e-9) * 1e3) / sec["dmax"] * NR, 0, NR - 1).astype(int)
    col_idx = np.clip(np.round(np.interp(xs, mdl.x, np.arange(mdl.nx))).astype(int), 0, mdl.nx - 1)
    dep_rgb = sec["img"][row, col_idx[None, :], :3]
    below = zeta_b < 0
    depth_px = np.where(below, -zeta_b / dz, 0)
    base = np.array([0.42, 0.44, 0.45]) * (1 - 0.3 * np.clip(depth_px / 60, 0, 1))[..., None]
    base = np.where((depth_px < 2.0)[..., None], np.array([0.78, 0.80, 0.81]), base)
    img = np.clip((bg * (1 - alpha[..., None]) + lit * alpha[..., None]) ** 0.95, 0, 1)
    img = np.where(in_dep[..., None], dep_rgb, img)
    img = np.where(below[..., None], base, img)
    if return_alpha:
        a_out = np.where(below | in_dep, 1.0, np.clip(alpha * 1.15, 0, 1))
        rgb = np.where((below | in_dep)[..., None], img, np.clip(lit ** 0.95, 0, 1))
        return np.dstack([rgb, a_out])
    return img


# ================================================================ Jupyter UI
try:
    import ipywidgets as W
    from IPython.display import display, clear_output
    _HAVE_WIDGETS = True
except ImportError:
    _HAVE_WIDGETS = False

C_SILI, C_CARB, C_VEL, C_FRONT, C_WATER = "#e0a526", "#1f6fb2", "#0d7f8c", "#d62728", "#3b8fc4"


def _fmt(v):
    a = abs(v)
    if v == 0 or not np.isfinite(v):
        return "0"
    if 0.01 <= a < 1000:
        return f"{v:.3g}"
    e = int(np.floor(np.log10(a)))
    return f"{v / 10**e:.2f}×10<sup>{e}</sup>"


class GeoFluxSimApp:
    VE3, DEPTH_X, DEPTH_Y = 3.0, -0.6, 0.30

    def __init__(self):
        self.model = TurbidityModel()
        self.ghosts, self.next_ghost = [], 2.0
        self.playing, self._task = False, None
        self._overlays = []
        be = matplotlib.get_backend().lower()
        self._ipympl = "ipympl" in be or "widget" in be
        self._build_widgets()
        self._build_figures()
        self.model.advance(12.0)
        self._redraw_static()
        self.render()

    # ----------------------------------------------------------- widgets
    def _build_widgets(self):
        st, lay = {"description_width": "150px"}, W.Layout(width="340px")
        def fs(v, a, b, s, d, fmt=".2f", width=None):
            return W.FloatSlider(value=v, min=a, max=b, step=s, description=d, style=st,
                                 layout=W.Layout(width=width) if width else lay,
                                 readout_format=fmt, continuous_update=False)
        def dd(opts, v, d):
            return W.Dropdown(options=list(opts), value=v, description=d, style=st, layout=lay)
        p = DEFAULTS
        self.w = dict(
            D_s=fs(p["D_s"], 10, 400, 5, "Median grain size (µm)", ".0f"),
            rho_s=fs(p["rho_s"], 2500, 2900, 10, "Grain density (kg/m³)", ".0f"),
            shape_s=dd(SH_SILI, p["shape_s"], "Grain shape"),
            D_c=fs(p["D_c"], 10, 400, 5, "Median grain size (µm)", ".0f"),
            rho_c=fs(p["rho_c"], 1500, 2850, 25, "Grain density (kg/m³)", ".0f"),
            shape_c=dd(SH_CARB, p["shape_c"], "Grain shape"),
            fc=fs(p["fc"], 0, 100, 5, "Carbonate fraction (%)", ".0f"),
            sort=fs(p["sort"], 1.0, 2.5, 0.1, "Sorting σg", ".1f"),
            C0=fs(p["C0"], 0.5, 15, 0.5, "Inflow conc. (vol %)", ".1f"),
            T_w=fs(p["T_w"], 0, 40, 1, "Water temperature (°C)", ".0f"),
            rho_a=fs(p["rho_a"], 995, 1030, 1, "Ambient density (kg/m³)", ".0f"),
            slope=fs(p["slope"], 0, 15, 0.5, "Bed slope (%)", ".1f"),
            U0=fs(p["U0"], 0.05, 0.5, 0.01, "Inflow velocity (m/s)"),
            h0=fs(p["h0"], 0.02, 0.20, 0.01, "Inflow thickness (m)"),
            T_in=fs(p["T_in"], 2, 60, 1, "Inflow duration (s)", ".0f"),
            Cd=fs(p["Cd"], 0.001, 0.02, 0.0005, "Bed drag Cd", ".4f"),
            ea=fs(p["ea"], 0.0, 0.15, 0.005, "Entrainment coeff. a", ".3f"),
            r0=fs(p["r0"], 1.0, 3.0, 0.1, "Near-bed ratio r0", ".1f"),
        )
        self.nu_lbl = W.HTML()
        # run bar
        self.play_btn = W.ToggleButton(value=False, description="▶ Play", button_style="success", layout=W.Layout(width="100px"))
        self.reset_btn = W.Button(description="⟲ Reset", button_style="warning", layout=W.Layout(width="100px"))
        self.probe_w = fs(3.0, 0.3, 7.7, 0.1, "Probe position (m)", ".1f", "330px")
        self.speed = fs(0.5, 0.1, 2.0, 0.1, "Sim s per frame", ".1f", "300px")
        self.status = W.HTML()
        # view options
        tb = {"description_width": "120px", "button_width": "105px"}
        self.flow_mode = W.ToggleButtons(options=["Realistic", "Concentration", "Velocity"], description="Show flow as", style=tb)
        self.mode_w = W.ToggleButtons(options=["Composition", "Grain size"], description="Deposit colours", style=tb)
        self.show_flow = W.Checkbox(value=True, description="Flow", indent=False, layout=W.Layout(width="70px"))
        self.show_dep = W.Checkbox(value=True, description="Deposit", indent=False, layout=W.Layout(width="80px"))
        self.show3d = W.Checkbox(value=True, description="Update 3D panel", indent=False, layout=W.Layout(width="140px"))
        self.live_plots = W.Checkbox(value=True, description="Update all plots while playing", indent=False, layout=W.Layout(width="230px"))
        self.exag_w = fs(5, 1, 40, 1, "Deposit exaggeration", ".0f", "330px")
        self.por_w = fs(0.40, 0.25, 0.6, 0.01, "Deposit porosity λ", ".2f", "330px")
        self.out_views, self.out_plots = W.Output(), W.Output()
        self.dep_stats, self.eq_html = W.HTML(), W.HTML()
        # events
        self.play_btn.observe(self._on_play, "value")
        self.reset_btn.on_click(lambda _: self.reset())
        for wid in self.w.values():
            wid.observe(lambda ch: self.reset(), "value")
        for wid in (self.probe_w, self.exag_w, self.por_w, self.mode_w, self.flow_mode,
                    self.show_flow, self.show_dep, self.show3d):
            wid.observe(self._on_display, "value")
        # layout: run bar, then PARAMETERS first, then views, then every plot
        hdr = lambda t: W.HTML(f"<b style='font-size:14px'>{t}</b>")
        box = lambda *c: W.VBox(list(c), layout=W.Layout(border="1px solid #d5e1e8", padding="6px 8px", margin="3px"))
        k = self.w
        params = W.GridBox([
            box(hdr("Siliciclastic sediment"), k["D_s"], k["rho_s"], k["shape_s"]),
            box(hdr("Carbonate sediment"), k["D_c"], k["rho_c"], k["shape_c"]),
            box(hdr("Mixture"), k["fc"], k["sort"], k["C0"]),
            box(hdr("Fluid"), k["T_w"], self.nu_lbl, k["rho_a"]),
            box(hdr("Flow and flume"), k["slope"], k["U0"], k["h0"], k["T_in"]),
            box(hdr("Model coefficients"), k["Cd"], k["ea"], k["r0"]),
        ], layout=W.Layout(grid_template_columns="repeat(auto-fill, minmax(360px, 1fr))"))
        header = W.HTML("<h3 style='margin:4px 0'>GeoFluxSim with David — turbidity currents, deposits and physics</h3>")
        runbar = W.HBox([self.play_btn, self.reset_btn, self.probe_w, self.speed, self.status],
                        layout=W.Layout(flex_flow="row wrap", align_items="center"))
        accordion = W.Accordion(children=[params])
        accordion.set_title(0, "Parameters (changing one restarts the run)")
        accordion.selected_index = 0
        viewbar = W.VBox([W.HBox([self.flow_mode, self.show_flow, self.show_dep, self.show3d],
                                 layout=W.Layout(flex_flow="row wrap", align_items="center")),
                          W.HBox([self.mode_w, self.exag_w, self.por_w, self.live_plots],
                                 layout=W.Layout(flex_flow="row wrap", align_items="center"))])
        self.ui = W.VBox([header, runbar, accordion, viewbar, self.out_views,
                          W.HTML("<h3 style='margin:10px 0 0'>Flow, deposits and physics at a glance</h3>"),
                          self.out_plots, self.dep_stats, self.eq_html])

    # ----------------------------------------------------------- figures
    def _build_figures(self):
        if not self._ipympl:
            plt.ioff()
        f = self.figA = plt.figure(figsize=(15, 9.4), dpi=75, constrained_layout=True)
        gs = f.add_gridspec(2, 3, height_ratios=[1, 1.35])
        self.ax_side = f.add_subplot(gs[0, :])
        self.ax_3d = f.add_subplot(gs[1, :2])
        self.ax_cam = f.add_subplot(gs[1, 2])
        f = self.figB = plt.figure(figsize=(15, 13.5), dpi=75, constrained_layout=True)
        gs = f.add_gridspec(3, 8, height_ratios=[1, 1.15, 1.25])
        self.ax_u = f.add_subplot(gs[0, 0:2]); self.ax_c = f.add_subplot(gs[0, 2:4])
        self.ax_x = f.add_subplot(gs[0, 4:6]); self.ax_bud = f.add_subplot(gs[0, 6:8])
        self.ax_sec = f.add_subplot(gs[1, 0:3])
        self.ax_core = f.add_subplot(gs[1, 3]); self.ax_core_d = f.add_subplot(gs[1, 4:6], sharey=self.ax_core)
        self.ax_dist = f.add_subplot(gs[1, 6:8])
        self.ax_core_t, self.ax_dist_t = self.ax_core_d.twiny(), self.ax_dist.twinx()   # created once, reused
        self.ax_ph = f.add_subplot(gs[2, 0:5]); self.ax_ws = f.add_subplot(gs[2, 5:8])
        if not self._ipympl:
            plt.close(self.figA); plt.close(self.figB)

    def _redraw_static(self):
        m = self.model
        self.ax_side.clear()
        self.side_ext = (0, m.L, -0.03, m.S * m.L + 0.32)
        self.im_side = self.ax_side.imshow(np.zeros((2, 2, 3)), extent=self.side_ext, aspect="auto", interpolation="bilinear")
        self.side_probe = self.ax_side.axvline(0, color=C_FRONT, lw=1, alpha=0.7)
        self.side_text = self.ax_side.text(0.99, 0.94, "", transform=self.ax_side.transAxes, ha="right", va="top",
                                           fontsize=11, bbox=dict(fc="white", ec="none", alpha=0.8))
        self.ax_side.set_xlabel("Distance downstream (m)"); self.ax_side.set_ylabel("Elevation (m)")
        self.cax = self.ax_side.inset_axes([0.012, 0.80, 0.18, 0.05])   # colour bar for field views
        self.ax_cam.clear()
        self.im_cam = self.ax_cam.imshow(np.zeros((2, 2, 3)), aspect="auto", interpolation="bilinear")
        self.ax_cam.set_xticks([]); self.ax_cam.set_yticks([])
        self.ax_cam.set_title("Probe camera (20 cm view)")

    # --------------------------------------------------------- helpers
    def _vmax(self):
        m = self.model
        return max(m.Cin.sum() * m.p["r0"], 1e-5), max(self.u_run, 1e-3) / U_NORM

    def _flow_kw(self):
        vc, vu = self._vmax()
        return dict(mode=self.flow_mode.value, show_flow=self.show_flow.value, show_dep=self.show_dep.value,
                    vmax_c=vc, vmax_u=vu)

    def _probe_profile(self, sec):
        """Height above (deposit-raised) bed and normalised field value at the probe."""
        m, xp = self.model, self.probe_w.value
        pr = m.probe(xp)
        if pr["h"] < 0.004 or pr["Ct"] < 1e-4 or not self.show_flow.value:
            return None
        vc, vu = self._vmax()
        zeta = np.linspace(0, 2.2 * pr["h"], 80)
        e = zeta / pr["h"]
        if self.flow_mode.value == "Concentration":
            val = np.clip(pr["Ct"] * m.ap * np.exp(-m.kp * e) * 0.5 * (1 - np.tanh((e - 1) / 0.12)) / vc, 0, 1)
        else:
            val = np.clip(pr["U"] * u_shape(e) / vu, 0, 1)
        zbase = m.S * (m.L - xp) + (np.interp(xp, m.x, sec["thick"]) * 1e-3 * self.exag_w.value if self.show_dep.value else 0)
        return zeta, val, zbase

    def _overlay_profile(self, ax, sec, xscale, zscale=1.0, lw=2.0):
        prof = self._probe_profile(sec)
        if prof is None or self.flow_mode.value == "Realistic":
            return
        zeta, val, zbase = prof
        xp = self.probe_w.value
        z = (zbase + zeta) * zscale
        self._overlays += ax.plot(xp + val * xscale, z, color="white", lw=lw, zorder=6)
        self._overlays += ax.plot([xp, xp], [z[0], z[-1]], color="white", lw=lw * 0.6, zorder=6)

    # ------------------------------------------------------------ views
    def _proj(self, x, y, z):
        d = (y + self.model.p["W"] / 2) / self.model.p["W"]
        return x + self.DEPTH_X * d, self.VE3 * z + self.DEPTH_Y * d

    def _draw_views(self, sec):
        m, ex, xp = self.model, self.exag_w.value, self.probe_w.value
        for a in self._overlays:
            a.remove()
        self._overlays = []
        kw = self._flow_kw()
        self.im_side.set_data(render_flow(m, sec, ex, *self.side_ext, 900, 210, **kw))
        mode = self.flow_mode.value
        self.ax_side.set_title(f"Flume side view — flow shown as {mode.lower()}"
                               + (f"; deposit ×{ex:.0f}" if self.show_dep.value else ""))
        self.side_probe.set_xdata([xp] * 2)
        self.side_text.set_text(f"t = {m.t:5.1f} s   front = {m.front:4.2f} m")
        self._overlay_profile(self.ax_side, sec, xscale=0.9)
        zb = m.S * (m.L - xp) + (np.interp(xp, m.x, sec["thick"]) * 1e-3 * ex if self.show_dep.value else 0)
        ext = (xp - 0.1, xp + 0.1, zb - 0.02, zb + 0.14)
        self.im_cam.set_data(render_flow(m, sec, ex, *ext, 320, 300, **kw))
        self.im_cam.set_extent(ext); self.ax_cam.set_xlim(ext[:2]); self.ax_cam.set_ylim(ext[2:])
        self._overlay_profile(self.ax_cam, sec, xscale=0.085, lw=2.5)
        # colour bar for field views
        self.cax.clear()
        if mode != "Realistic" and self.show_flow.value:
            vc, vu = self._vmax()
            vmax = vc * 100 if mode == "Concentration" else vu * 100
            grad = np.linspace(0, 1, 128)[None, :]
            self.cax.imshow(field_cmap(mode)(0.12 + 0.86 * grad if mode == "Concentration" else grad), aspect="auto",
                            extent=(0, vmax, 0, 1))
            self.cax.set_yticks([]); self.cax.set_xticks([0, vmax])
            self.cax.set_xticklabels(["0", f"{vmax:.1f} " + ("vol %" if mode == "Concentration" else "cm/s")], fontsize=9)
            self.cax.set_title(f"{mode} (white line: profile at probe)", fontsize=9, loc="left")
            self.cax.set_visible(True)
        else:
            self.cax.set_visible(False)
        if self.show3d.value or not self.ax_3d.collections:
            self._draw_3d(sec)

    def _draw_3d(self, sec):
        m, ax, Wd, ex = self.model, self.ax_3d, self.model.p["W"], self.exag_w.value
        ax.clear()
        ztop = m.S * m.L + 0.30
        ax.imshow(np.linspace(1, 0, 64)[:, None, None] * WATER_TOP + np.linspace(0, 1, 64)[:, None, None] * WATER_DEEP,
                  extent=(-0.8, m.L + 0.2, -0.15, self.VE3 * ztop + 0.55), aspect="auto", zorder=0)
        for (x0, z0), (x1, z1) in [((0, m.S * m.L), (m.L, 0)), ((0, ztop), (m.L, 0.30)), ((0, m.S * m.L), (0, ztop)), ((m.L, 0), (m.L, 0.30))]:
            px, py = self._proj(np.array([x0, x1]), Wd / 2, np.array([z0, z1]))
            ax.plot(px, py, color="#e9f4fa", lw=1, alpha=0.8, zorder=1)
        xg, yg = np.linspace(0, m.L, 160), np.linspace(-Wd / 2, Wd / 2, 12)
        X, Y = np.meshgrid(xg, yg)
        h = np.interp(xg, m.x, m.h)
        C = np.interp(xg, m.x, m.m.sum(0) / m.h)
        on = (C > 1e-3 * m.Cin.sum()) & (h > 0.004)
        n = _bilinear(NOISE2, ((X - m.phase) / 4.0 * 512) % 512, (Y / Wd * 40 + m.t * 0.6) % 64)
        dh = np.interp(xg, m.x, sec["thick"]) * 1e-3 * ex * self.show_dep.value
        on = on & self.show_flow.value
        Z = m.S * (m.L - X) + dh[None, :] + np.where(on[None, :], h[None, :] * (0.8 + 0.7 * (n - 0.5)), 0.0)
        nrm = np.dstack([-np.gradient(Z * self.VE3, xg, axis=1), -np.gradient(Z * self.VE3, yg, axis=0), np.ones_like(Z)])
        nrm /= np.linalg.norm(nrm, axis=2, keepdims=True)
        lamb = np.clip(nrm @ np.array([0.35, 0.30, 0.89]), 0, 1)
        top = np.array([np.interp(xg, m.x, sec["top"][:, j]) for j in range(3)]).T
        bed = np.where((dh > 1e-5)[:, None], top, np.array([0.36, 0.37, 0.37]))
        cfrac = np.where(on, np.clip(C / max(m.Cin.sum(), 1e-9) * 1.5, 0.35, 1), 0)
        cur = np.array([0.80, 0.73, 0.60]) * (0.75 + 0.5 * n[..., None])
        if self.flow_mode.value != "Realistic":
            vc, vu = self._vmax()
            Ug = np.interp(xg, m.x, m.vel(m.h, m.q))
            fv = np.clip(C * m.p["r0"] / vc, 0, 1) if self.flow_mode.value == "Concentration" else np.clip(Ug / U_NORM / vu, 0, 1)
            cm = field_cmap(self.flow_mode.value)(0.12 + 0.86 * fv if self.flow_mode.value == "Concentration" else fv)[:, :3]
            cur = np.broadcast_to(cm[None, :, :], X.shape + (3,)) * (0.85 + 0.3 * n[..., None])
        base = bed[None] * (1 - cfrac)[None, :, None] + cur * cfrac[None, :, None]
        base = base * (1 - 0.25 * ((Y + Wd / 2) / Wd))[..., None]
        col = np.clip(base * (0.30 + 0.85 * lamb[..., None]), 0, 1)
        PX, PY = self._proj(X, Y, Z)
        quads = np.stack([np.dstack([PX[:-1, :-1], PY[:-1, :-1]]), np.dstack([PX[:-1, 1:], PY[:-1, 1:]]),
                          np.dstack([PX[1:, 1:], PY[1:, 1:]]), np.dstack([PX[1:, :-1], PY[1:, :-1]])], axis=2)
        fc = 0.25 * (col[:-1, :-1] + col[:-1, 1:] + col[1:, 1:] + col[1:, :-1])
        ax.add_collection(PolyCollection(quads[::-1].reshape(-1, 4, 2), facecolors=fc[::-1].reshape(-1, 3),
                                         edgecolors=fc[::-1].reshape(-1, 3), linewidths=0.3, zorder=2))
        rgba = render_flow(m, sec, ex, 0, m.L, -0.02, ztop, 700, 170, return_alpha=True, **self._flow_kw())
        ax.imshow(rgba, extent=(0, m.L, self.VE3 * -0.02, self.VE3 * ztop), aspect="auto", interpolation="bilinear", zorder=3)
        for (x0, z0), (x1, z1) in [((0, ztop), (m.L, 0.30)), ((0, -0.02), (m.L, -0.02)), ((0, -0.02), (0, ztop)), ((m.L, -0.02), (m.L, 0.30))]:
            ax.plot([x0, x1], [self.VE3 * z0, self.VE3 * z1], color="#f2f9fc", lw=1.4, alpha=0.9, zorder=4)
        for xc, zc in [(0, ztop), (m.L, 0.30)]:
            px, py = self._proj(np.array([xc, xc]), np.array([-Wd / 2, Wd / 2]), np.array([zc, zc]))
            ax.plot(px, py, color="#f2f9fc", lw=1.2, alpha=0.9, zorder=4)
        xp = self.probe_w.value
        zp = m.S * (m.L - xp) + np.interp(xp, m.x, sec["thick"]) * 1e-3 * ex * self.show_dep.value
        px, py = self._proj(np.array([xp, xp]), np.array([-Wd / 2, Wd / 2]), np.array([zp, zp]))
        ax.plot(px, py, color=C_FRONT, lw=2, zorder=2.5)
        ax.set_xlim(-0.8, m.L + 0.2); ax.set_ylim(-0.15, self.VE3 * ztop + 0.55); ax.set_axis_off()
        prof = self._probe_profile(sec)
        if prof is not None and self.flow_mode.value != "Realistic":
            zeta, val, zbase = prof
            ax.plot(xp + val * 0.9, self.VE3 * (zbase + zeta), color="white", lw=2, zorder=6)
            ax.plot([xp, xp], [self.VE3 * zbase, self.VE3 * (zbase + zeta[-1])], color="white", lw=1.2, zorder=6)
        ax.set_title("3D view through the near glass wall (vertical ×3; deposit layers show on the wall)")

    # ------------------------------------------------------------ plots
    def _draw_flow(self):
        m, xp = self.model, self.probe_w.value
        pr = m.probe(xp)
        ztop = max(2.4 * m.p["h0"], 2.0 * pr["h"], 0.05) * 100
        z = np.linspace(0, ztop, 200)
        for ax in (self.ax_u, self.ax_c, self.ax_x):
            ax.clear(); ax.grid(alpha=0.3)
        has = pr["h"] > 0.003 and pr["Ct"] > 1e-5
        for k, gst in enumerate(self.ghosts):
            if gst["h"] > 0.003:
                self.ax_u.plot(gst["U"] * u_shape(z / 100 / gst["h"]) * 100, z, color=C_VEL, alpha=0.15 + 0.07 * k, lw=1)
        if has:
            self.ax_u.plot(pr["U"] * u_shape(z / 100 / pr["h"]) * 100, z, color=C_VEL, lw=2.4)
            self.ax_u.axhline(ETA_M * pr["h"] * 100, color="grey", ls=":", lw=1)
            e = z / 100 / pr["h"]
            cs = m.ap * np.exp(-m.kp * e) * 0.5 * (1 - np.tanh((e - 1) / 0.12))
            self.ax_c.plot(pr["Ct"] * cs * 100, z, color="k", lw=2.4, label="Total")
            self.ax_c.plot(pr["Cs"] * cs * 100, z, color=C_SILI, lw=2, label="Siliciclastic")
            self.ax_c.plot(pr["Cc"] * cs * 100, z, color=C_CARB, lw=2, label="Carbonate")
            self.ax_c.legend(loc="upper right", fontsize=9)
        else:
            for ax in (self.ax_u, self.ax_c):
                ax.text(0.5, 0.5, "No current at the probe yet", transform=ax.transAxes, ha="center", color="grey")
        umax = 1.35 * max([m.p["U0"], abs(pr["U"])] + [abs(gst["U"]) for gst in self.ghosts]) / U_NORM * 100
        self.ax_u.set_xlim(0, max(umax, 5)); self.ax_u.set_ylim(0, ztop)
        self.ax_c.set_xlim(0, max(m.Cin.sum() * m.p["r0"], pr["Ct"] * m.p["r0"], 1e-3) * 105); self.ax_c.set_ylim(0, ztop)
        self.ax_u.set(xlabel="Velocity (cm/s)", ylabel="Height above bed (cm)", title="Velocity profile at the probe")
        self.ax_c.set(xlabel="Concentration (vol %)", ylabel="Height above bed (cm)", title="Concentration profile at the probe")
        self.ax_x.plot(m.t_hist, m.xf_hist, color=C_FRONT, lw=2)
        self.ax_x.axhline(xp, color=C_FRONT, ls="--", lw=1, alpha=0.4)
        self.ax_x.set(xlim=(0, m.p["t_max"]), ylim=(0, m.L * 1.02), xlabel="Time (s)", ylabel="Front position (m)",
                      title="Front position over time")
        ax = self.ax_bud; ax.clear()
        b = np.array(m.budget)
        ax.plot(b[:, 0], b[:, 1], color=C_WATER, lw=2, label="In suspension")
        ax.plot(b[:, 0], b[:, 2], color=C_SILI, lw=2, label="Deposited")
        ax.plot(b[:, 0], b[:, 3], color=C_FRONT, lw=2, label="Bypassed outlet")
        ax.set(xlim=(0, m.p["t_max"]), ylim=(0, 100), xlabel="Time (s)", ylabel="Share of injected sediment (%)",
               title="Sediment budget over time")
        ax.legend(fontsize=9); ax.grid(alpha=0.3)


    def _draw_deposits(self, sec):
        m, xp = self.model, self.probe_w.value
        for ax in (self.ax_sec, self.ax_core, self.ax_core_d, self.ax_dist, self.ax_core_t, self.ax_dist_t):
            ax.clear()
        ax = self.ax_sec
        ymax = sec["dmax"] * 1.3
        if sec["layers"] is not None and sec["thick"].max() > 0:
            ax.imshow(sec["img"], origin="lower", extent=(0, m.L, 0, sec["dmax"]), aspect="auto", interpolation="nearest")
        else:
            ax.text(0.5, 0.5, "No deposit yet. Press Play.", transform=ax.transAxes, ha="center", color="grey")
        ax.axvline(xp, color=C_FRONT, lw=1.5)
        r = sec["rng"]
        span = (f"carbonate {r[0] * 100:.0f}% (amber) → {r[1] * 100:.0f}% (blue)" if self.mode_w.value == "Composition"
                else f"mean size {np.exp(r[2]):.0f} µm (pale) → {np.exp(r[3]):.0f} µm (dark)")
        ax.set(xlim=(0, m.L), ylim=(0, ymax), xlabel="Distance downstream (m)", ylabel="Deposit thickness (mm)",
               title="Deposit cross-section along the flume")
        ax.text(0.99, 0.97, f"colours stretched to this deposit: {span}\ndark lines: time lines every 7.5 s",
                transform=ax.transAxes, ha="right", va="top", fontsize=9, bbox=dict(fc="white", ec="none", alpha=0.8))
        ax.grid(alpha=0.25)
        # core
        i = int(np.clip(round(xp / m.dx - 0.5), 0, m.nx - 1))
        T = sec["thick"][i]
        if sec["layers"] is not None and T > 1e-4:
            th, cf, Dm = (a[:, i] for a in sec["layers"])
            keep = th > 0
            th, cf, Dm = th[keep], cf[keep], Dm[keep]
            zc = np.cumsum(th); zmid = zc - th / 2
            cols = deposit_colour(cf, np.log(Dm), sec["rng"], self.mode_w.value)
            self.ax_core.bar(np.zeros_like(th), th, bottom=zc - th, width=0.6, color=cols, edgecolor="none")
            self.ax_core.set_xlim(-0.5, 0.5); self.ax_core.set_xticks([])
            self.ax_core_d.plot(Dm, zmid, color=C_VEL, lw=2, label="Mean grain size (µm)")
            tw = self.ax_core_t
            tw.xaxis.tick_top(); tw.xaxis.set_label_position("top")
            tw.plot(cf * 100, zmid, color=C_CARB, lw=2)
            tw.set_xlim(0, 100); tw.set_xlabel("Carbonate (%)", color=C_CARB)
            self.ax_core_d.set_xlabel("Mean grain size (µm)", color=C_VEL)
            self.ax_core_d.set_xlim(np.exp(sec["lnDmin"]) * 0.9, np.exp(sec["lnDmax"]) * 1.1)
        else:
            self.ax_core.text(0.5, 0.5, "No deposit\nat the probe", transform=self.ax_core.transAxes, ha="center", color="grey")
        self.ax_core.set(ylim=(0, max(T * 1.15, 0.05)), ylabel="Height in deposit (mm)", title=f"Core at {xp:.1f} m")
        self.ax_core_d.set_title("Grading in the core"); self.ax_core_d.grid(alpha=0.3)
        # downstream distribution
        ax = self.ax_dist
        ax.fill_between(m.x, 0, sec["mmS"], color=C_SILI, lw=0, label="Siliciclastic")
        ax.fill_between(m.x, sec["mmS"], sec["thick"], color=C_CARB, lw=0, label="Carbonate")
        ax.axvline(xp, color=C_FRONT, lw=1.2)
        ax.set(xlim=(0, m.L), ylim=(0, ymax), xlabel="Distance downstream (m)", ylabel="Thickness (mm)",
               title="Downstream distribution")
        ax.legend(loc="upper right", fontsize=9); ax.grid(alpha=0.3)
        tw = self.ax_dist_t
        tw.yaxis.tick_right(); tw.yaxis.set_label_position("right")
        good = sec["thick"] > 0.005 * sec["dmax"]
        tw.plot(m.x[good], sec["meanD"][good], color=C_VEL, lw=2)
        tw.set_ylim(np.exp(sec["lnDmin"]) * 0.9, np.exp(sec["lnDmax"]) * 1.1)
        tw.set_ylabel("Mean grain size (µm)", color=C_VEL)
        # stats
        thick, tot = sec["thick"], sec["thick"].sum()
        if tot > 0:
            thirds = [thick[j * m.nx // 3:(j + 1) * m.nx // 3].sum() / tot * 100 for j in range(3)]
            imax = int(np.argmax(thick))
            run = m.x[np.nonzero(thick > 0.01 * sec["dmax"])[0][-1]]
            self.dep_stats.value = (
                "<div style='display:flex;gap:24px;flex-wrap:wrap;font-family:sans-serif;font-size:14px;padding:6px 4px'>"
                f"<div>Thickest point<br><b>{thick[imax]:.2f} mm at {m.x[imax]:.1f} m</b></div>"
                f"<div>Deposit centre of mass<br><b>{(thick * m.x).sum() / tot:.2f} m</b></div>"
                f"<div>Runout (1% of max)<br><b>{run:.2f} m</b></div>"
                f"<div>Proximal / medial / distal thirds<br><b>{thirds[0]:.0f}% / {thirds[1]:.0f}% / {thirds[2]:.0f}%</b></div>"
                f"<div>Carbonate share of deposit<br><b>{sec['mmC'].sum() / tot * 100:.0f}%</b></div></div>")
        else:
            self.dep_stats.value = ""

    def _draw_physics(self):
        m, ph = self.model, self.model.physics(self.probe_w.value)
        ax = self.ax_ph
        ax.clear(); ax.set_xlim(0, 520); ax.set_ylim(300, 0); ax.set_axis_off()
        cx, cy, a = 260, 170, np.radians(7)
        def R(x, y):
            x, y = np.asarray(x, float), np.asarray(y, float)
            return cx + (x - cx) * np.cos(a) - (y - cy) * np.sin(a), cy + (x - cx) * np.sin(a) + (y - cy) * np.cos(a)
        def poly(xs, ys, **kw):
            px, py = R(xs, ys); ax.fill(px, py, **kw)
        def arrow(x1, y1, x2, y2, color, lw=3):
            (a1, b1), (a2, b2) = R(x1, y1), R(x2, y2)
            ax.annotate("", xy=(a2, b2), xytext=(a1, b1),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=12 + lw * 2))
        def text(x, y, s, **kw):
            px, py = R(x, y); ax.text(px, py, s, rotation=-7, rotation_mode="anchor", fontsize=11, **kw)
        ax.imshow(np.linspace(0, 1, 32)[:, None, None] * np.array([0.50, 0.71, 0.84]) +
                  np.linspace(1, 0, 32)[:, None, None] * np.array([0.78, 0.90, 0.96]), extent=(0, 520, 300, 0), aspect="auto")
        bedY, depY, topY = 200, 188, 100
        poly([-60, 600, 600, -60], [bedY, bedY, 400, 400], color="#5d6264")
        poly([-60, 600, 600, -60], [depY, depY, bedY, bedY], color=C_SILI)
        poly([-60, 600, 600, -60], [depY, depY, depY + 5, depY + 5], color=C_CARB)
        xs = np.linspace(-60, 600, 80)
        poly(np.r_[xs, 600, -60], np.r_[topY + 5 * np.sin(xs / 9) + 3 * np.sin(xs / 4), depY, depY], color="#dacdb3")
        if ph["active"]:
            for j in range(1, 9):
                e = j / 8 * 1.3; zz = depY - e / 1.3 * (depY - topY - 4); u = float(u_shape(e) / u_shape(ETA_M))
                if u > 0.04:
                    arrow(20, zz, 20 + u * 60, zz, C_VEL, 1.4)
            text(16, topY - 12, "u(z)", color="#14232c")
            mom = [("Gravity drive", ph["grav"], "#c0392b"), ("Pressure gradient", ph["press"], "#8e44ad"),
                   ("Entrainment drag", ph["entr"], C_WATER), ("Bed friction", ph["fric"], "#1d2b33")]
            mmax = max(max(abs(v) for _, v, _ in mom), 1e-12)
            for kk, (name, v, colr) in enumerate(mom):
                ln = 16 + 150 * abs(v) / mmax; y = [122, 142, 162, 180][kk]; x0 = 250 - ln / 2
                arrow(x0, y, x0 + ln, y, colr) if v >= 0 else arrow(x0 + ln, y, x0, y, colr)
                text(x0 + ln + 8, y + 4, name, color="#14232c")
            ew_r = ph["ew"] * abs(ph["U"])
            for x in (420, 455, 490):
                arrow(x, topY - 48, x, topY - 6, C_WATER, min(1.5 + ew_r * 2500, 6))
            text(400, topY - 54, "Water entrained", color="#14232c")
            smax = max(ph["depS"], ph["depC"], 1e-12)
            arrow(104, topY + 24, 104, depY - 2, C_SILI, 1.5 + 6 * ph["depS"] / smax)
            arrow(134, topY + 24, 134, depY - 2, C_CARB, 1.5 + 6 * ph["depC"] / smax)
            text(96, topY + 16, "Settling", color="#14232c")
            ax.text(14, 24, f"Downstream →   Fr {ph['Fr']:.3g} ({'supercritical' if ph['Fr'] > 1 else 'subcritical'})   "
                            f"Ri {ph['Ri']:.3g}   Re {ph['Re']:,.0f}", fontsize=12, color="#14232c")
        else:
            ax.text(260, 40, "No current at the probe: move the probe or press Play", ha="center", fontsize=13, color="#14232c")
        ax.set_title("Forces and fluxes at the probe (arrow lengths scale with the live values)")
        # settling curves
        ax = self.ax_ws; ax.clear()
        Dg = np.logspace(1, 3, 120); p = m.p
        ax.loglog(Dg, stokes_velocity(Dg, p["rho_s"], m.nu, p["rho_a"]) * 1e3, "--", color="grey", label="Stokes (spheres)")
        ax.loglog(Dg, settling_velocity(Dg, p["rho_s"], *SH_SILI[p["shape_s"]], m.nu, p["rho_a"]) * 1e3, color=C_SILI, lw=2, label="Siliciclastic")
        ax.loglog(Dg, settling_velocity(Dg, p["rho_c"], *SH_CARB[p["shape_c"]], m.nu, p["rho_a"]) * 1e3, color=C_CARB, lw=2, label="Carbonate")
        ax.scatter(m.D, np.maximum(m.ws0 * 1e3, 1e-3), c=[C_SILI if c == 0 else C_CARB for c in m.comp], edgecolor="white", zorder=3, s=45)
        ax.set(xlabel="Grain diameter (µm)", ylabel="Settling velocity (mm/s)", ylim=(1e-3, 1e3),
               title="Settling velocity vs grain size (dots: this run's six classes)")
        ax.legend(fontsize=9); ax.grid(alpha=0.3, which="both")
        # equation table
        por = self.por_w.value
        act = ph["active"]
        rows = [
            ("Excess density", "Δρ/ρ = Σ R<sub>i</sub>C<sub>i</sub>, R<sub>i</sub> = (ρ<sub>s,i</sub> − ρ<sub>w</sub>)/ρ<sub>w</sub>",
             f"{ph['B'] * 100:.3f}%", f"g′ = {ph['gp']:.3g} m/s²"),
            ("Gravity drive", "g′ h S", f"{_fmt(ph['grav'])} m²/s²", f"h = {ph['h'] * 100:.1f} cm, S = {p['slope']}%"),
            ("Pressure gradient", "−∂(½ g′h²)/∂x", f"{_fmt(ph['press'])} m²/s²", "pushes downstream" if ph["press"] >= 0 else "holds the flow back"),
            ("Bed friction", "−C<sub>d</sub> U|U|", f"{_fmt(ph['fric'])} m²/s²", f"U = {ph['U'] * 100:.1f} cm/s"),
            ("Water entrainment", "e<sub>w</sub> = a/√(1 + 718 Ri<sup>2.4</sup>), Ri = g′h/U²",
             f"{ph['ew'] * abs(ph['U']) * 1e3:.2f} mm/s", f"e<sub>w</sub> = {ph['ew']:.3g}; slows layer by {_fmt(ph['entr'])} m²/s²"),
            ("Densimetric Froude", "Fr = U/√(g′h)", f"{ph['Fr']:.3g}", "supercritical" if ph["Fr"] > 1 else "subcritical"),
            ("Settling (median sizes)", "w<sub>s</sub> = R g D² / (C<sub>1</sub>ν + √(0.75 C<sub>2</sub> R g D³))",
             f"{m.ws0[1] * 1e3:.2f} / {m.ws0[4] * 1e3:.2f} mm/s", f"siliciclastic / carbonate; ν = {m.nu * 1e6:.2f} mm²/s"),
            ("Hindered settling", "× (1 − C)<sup>4.65</sup>", f"× {ph['hind']:.3f}", f"C = {ph['pr']['Ct'] * 100:.2f} vol %"),
            ("Deposition rate", "Σ r<sub>0</sub> w<sub>s,i</sub> C<sub>i</sub> / (1 − λ)",
             f"{(ph['depS'] + ph['depC']) / (1 - por) * 6e4:.3f} mm/min",
             f"siliciclastic {ph['depS'] / (1 - por) * 6e4:.3f}, carbonate {ph['depC'] / (1 - por) * 6e4:.3f}"),
            ("Rouse number", "P = w<sub>s</sub>/(κu<sub>*</sub>), u<sub>*</sub> = √C<sub>d</sub>|U|",
             f"{ph['rouse'][0]:.3g} / {ph['rouse'][1]:.3g}",
             "; ".join(("siliciclastic " if k == 0 else "carbonate ") +
                       ("fully suspended" if r < 0.8 else "suspended, settling out" if r < 2.5 else "near-bed, deposits fast")
                       for k, r in enumerate(ph["rouse"]))),
            ("Flow Reynolds number", "Re = U h / ν", f"{ph['Re']:,.0f}", "turbulent" if ph["Re"] > 2000 else "laminar to transitional"),
        ]
        tr = "".join(f"<tr><td style='padding:4px 10px'><b>{a}</b></td><td style='padding:4px 10px'><i>{b}</i></td>"
                     f"<td style='padding:4px 10px;text-align:right'><b>{c if act or a.startswith('Settling') else '–'}</b></td>"
                     f"<td style='padding:4px 10px;color:#666'>{d if act or a.startswith('Settling') else ''}</td></tr>"
                     for a, b, c, d in rows)
        self.eq_html.value = ("<h4 style='margin:6px 0'>How each equation is applied at the probe</h4>"
                              f"<table style='font-family:sans-serif;font-size:13px;border-collapse:collapse'>{tr}</table>")

    # -------------------------------------------------------------- render
    def _show(self, fig, out):
        if self._ipympl:
            fig.canvas.draw_idle()
        else:
            with out:
                clear_output(wait=True)
                display(fig)

    def render(self, plots=True):
        m = self.model
        self.u_run = max(getattr(self, "u_run", m.p["U0"]), m.p["U0"])
        inside = m.m.sum(0) / m.h > 1e-3 * m.Cin.sum()
        if inside.any():
            self.u_run = max(self.u_run, float(np.abs(m.vel(m.h, m.q))[inside].max()))
        sec = build_section(m, self.por_w.value, self.mode_w.value)
        self._draw_views(sec)
        self._show(self.figA, self.out_views)
        if plots:
            self._draw_flow(); self._draw_deposits(sec); self._draw_physics()
            self._show(self.figB, self.out_plots)
        inj = max(m.injected, 1e-12)
        pr = m.probe(self.probe_w.value)
        self.nu_lbl.value = f"<span style='color:#666;font-size:12px;margin-left:155px'>ν = {m.nu * 1e6:.2f} mm²/s</span>"
        self.status.value = (
            f"<div style='padding:4px 12px;font-family:sans-serif;font-size:13px'>"
            f"<b>t</b> {m.t:.1f} s &nbsp; <b>front</b> {m.front:.2f} m &nbsp; <b>speed</b> {m.front_speed() * 100:.1f} cm/s &nbsp; "
            f"<b>probe h</b> {pr['h'] * 100 if pr['Ct'] > 1e-4 else 0:.1f} cm &nbsp; "
            f"<b>deposited</b> {100 * m.dep.sum() * m.dx / inj:.0f}% &nbsp; <b>bypassed</b> {100 * m.exited / inj:.0f}%"
            f"{' &nbsp; <b>— run complete</b>' if m.done else ''}</div>")

    # ------------------------------------------------------------ controls
    def reset(self):
        params = {k: w.value for k, w in self.w.items()}
        self.model = TurbidityModel(**params)
        self.ghosts, self.next_ghost = [], 2.0
        self.u_run = self.model.p["U0"]
        self._redraw_static()
        self.render()

    def _on_display(self, ch):
        if ch["owner"] is self.probe_w:
            self.ghosts = []
        self.render()

    def _on_play(self, change):
        self.playing = change["new"]
        self.play_btn.description = "⏸ Pause" if self.playing else "▶ Play"
        self.play_btn.button_style = "info" if self.playing else "success"
        if self.playing and (self._task is None or self._task.done()):
            self._task = asyncio.ensure_future(self._loop())
        if not self.playing:
            self.render()

    async def _loop(self):
        while self.playing:
            if self.model.done:
                self.play_btn.value = False
                break
            self.model.advance(self.speed.value)
            if self.model.t >= self.next_ghost:
                self.ghosts = (self.ghosts + [self.model.probe(self.probe_w.value)])[-6:]
                self.next_ghost += 2.0
            self.render(plots=self.live_plots.value)
            await asyncio.sleep(0.02)

    def show(self):
        display(self.ui)
        if self._ipympl:
            with self.out_views:
                display(self.figA.canvas)
            with self.out_plots:
                display(self.figB.canvas)
        return self


if _HAVE_WIDGETS:
    try:
        get_ipython  # noqa: F821
        app = GeoFluxSimApp().show()
    except NameError:
        print("Run this inside Jupyter for the interactive version.")
else:
    print("ipywidgets not found: run  pip install ipywidgets  and re-run this cell.")
