import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import norm
import cvxpy as cp

# ---------- RBF Histogram Matching via QP ----------
def match_histogram_rbf_qp(data, y_exp, sigma, num_centers=50, sigma_rbf=0.05):
    N = len(data)
    data = np.asarray(data)

    # Define RBF centers and basis matrix
    centers = np.linspace(np.min(data), np.max(data), num_centers)
    basis = np.exp(-0.5 * ((data[:, None] - centers[None, :]) / sigma_rbf)**2)

    # Target Gaussian histogram at RBF centers
    p_target = norm.pdf(centers, loc=y_exp, scale=sigma)
    p_target /= np.sum(p_target)

    # Define QP problem
    w = cp.Variable(N)
    objective = cp.Minimize(0.5 * cp.sum_squares(basis.T @ w - p_target))
    constraints = [w >= 0, cp.sum(w) == 1]
    prob = cp.Problem(objective, constraints)

    # Solve using OSQP with tight tolerances
    try:
        prob.solve(solver=cp.OSQP, eps_abs=1e-9, eps_rel=1e-9, max_iter=10_000)
    except cp.SolverError:
        print("⚠️ OSQP failed, trying ECOS as fallback...")
        prob.solve(solver=cp.ECOS, abstol=1e-8, reltol=1e-8, max_iters=10_000)

    if w.value is None:
        raise RuntimeError("QP solver failed to find a solution.")

    weights = np.clip(np.array(w.value).flatten(), 0, None)
    weights /= np.sum(weights)
    
    # Compute reweighted distribution
    p_reweighted = basis.T @ weights
    return weights, centers, p_target, p_reweighted

# ---------- CLI Wrapper ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile", help="Text file with y_exp, sigma, then simulated values")
    ap.add_argument("--out_prefix", default=None, help="Output file prefix")
    ap.add_argument("--num_centers", type=int, default=20)
    ap.add_argument("--sigma_rbf", type=float, default=0.01)
    args = ap.parse_args()

    with open(args.infile) as f:
        lines = [line.strip() for line in f if line.strip()]
    y_exp, sigma = map(float, lines[0].replace(",", " ").split())
    data = np.array([float(x) for x in lines[1:]])

    weights, centers, p_exp, p_w = match_histogram_rbf_qp(
        data, y_exp, sigma, num_centers=args.num_centers, sigma_rbf=args.sigma_rbf
    )

    out_prefix = args.out_prefix or Path(args.infile).with_suffix("").name
    np.savetxt(f"{out_prefix}_weights_rbf.txt", weights)
    np.save(f"{out_prefix}_weights_rbf.npy", weights)

    # Plotting
    plt.figure(figsize=(7, 4.5))
    #plt.hist(data, bins=50, density=True, histtype='step', label="Unweighted")
    plt.plot(centers, p_exp, label="Target Gaussian", linewidth=2)
    plt.step(centers, p_w, where='mid', linestyle='--', label=f"Reweighted (RBF)")
    plt.xlabel("Observable")
    plt.ylabel("Density")
    plt.title("Smoothed Histogram Matching via RBF + QP")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_rbf_match.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
