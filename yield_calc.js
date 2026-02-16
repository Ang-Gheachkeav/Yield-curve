// ==============================================================================
//  Yield Calculator for NBC
// ------------------------------------------------------------------------------
//  Version 	: 0.9.0
//  Created		: 2026-01-20
//  Author		: NBC MMD & JICA
//  Instrument	: Fixed-Rate Bond
// ============================================================================== 

//Find Left and Right Values from Array
function getBetween(days, xDays, yArr) {
	// Find Higest left including the same value
	let left = null;
	for (let i = 0; i < xDays.length; i++) {
		const d = xDays[i];
		const y = yArr[i];

		if (y == null) continue;
		if (d > days) continue; //Excluding the same value

		if (left === null || d >= left.day) { //Invluding the same value
			left = { idx: i, day: d, yield: y };
		}
	}

	// find Lowest right
	let right = null;
	for (let i = 0; i < xDays.length; i++) {
		const d = xDays[i];
		const y = yArr[i];
		if (y == null) continue;
		if (d <= days) continue;

		if (right === null || d < right.day) {
			right = { idx: i, day: d, yield: y };
		}
	}

	return { left, right };
}

class NSSFit {
	// ----- NSS basis functions (given t in years, tau1,tau2) -----
	static _factor(t, tau) {
		const x = t / tau;
		if (x === 0) return 1;
		return (1 - Math.exp(-x)) / x;
	}

	// Return basis vector b = [1, f1, f2, f3] so that y = beta · b
	static basis(t, tau1, tau2) {
		const f1 = NSSFit._factor(t, tau1);
		const e1 = Math.exp(-t / tau1);
		const f2 = f1 - e1;

		const f1b = NSSFit._factor(t, tau2);
		const e2 = Math.exp(-t / tau2);
		const f3 = f1b - e2;

		return [1, f1, f2, f3];
	}

	// ----- Utilities: build (t,y) samples skipping null/undefined -----
	static buildSamples(YearsArray, YieldArray) {
	const xs = [];
	const ys = [];
	for (let i = 0; i < YearsArray.length; i++) {
		const y = YieldArray[i];
		if (y == null || !Number.isFinite(y)) continue;
		const t = YearsArray[i];
		if (!Number.isFinite(t) || t <= 0) continue;
		xs.push(t);
		ys.push(y);
	}
	return { xs, ys };
	}

	// ----- Solve linear least squares for beta: minimize ||A*beta - y||^2 -----
	// A is n x 4, we solve normal equation: (A^T A) beta = A^T y
	static solveBetaLeastSquares(A, y) {
		// Build AtA (4x4) and AtY (4)
		const AtA = Array.from({ length: 4 }, () => Array(4).fill(0));
		const AtY = Array(4).fill(0);

		for (let i = 0; i < A.length; i++) {
			const r = A[i];
			for (let p = 0; p < 4; p++) {
			AtY[p] += r[p] * y[i];
			for (let q = 0; q < 4; q++) {
				AtA[p][q] += r[p] * r[q];
			}
			}
		}

		// Solve 4x4 linear system with Gaussian elimination (partial pivot)
		const M = AtA.map(row => row.slice());
		const b = AtY.slice();

		for (let col = 0; col < 4; col++) {
			// pivot
			let piv = col;
			for (let r = col + 1; r < 4; r++) {
				if (Math.abs(M[r][col]) > Math.abs(M[piv][col])) piv = r;
			}
			if (Math.abs(M[piv][col]) < 1e-12) return { ok: false, beta: null };

			// swap
			if (piv !== col) {
				[M[piv], M[col]] = [M[col], M[piv]];
				[b[piv], b[col]] = [b[col], b[piv]];
			}

			// eliminate
			for (let r = col + 1; r < 4; r++) {
			const f = M[r][col] / M[col][col];
				for (let c = col; c < 4; c++) M[r][c] -= f * M[col][c];
				b[r] -= f * b[col];
			}
		}

		// back-substitution
		const beta = Array(4).fill(0);
		for (let r = 3; r >= 0; r--) {
			let s = b[r];
			for (let c = r + 1; c < 4; c++) s -= M[r][c] * beta[c];
			beta[r] = s / M[r][r];
		}

		return { ok: true, beta };
	}

	// Compute SSE for given beta,tau1,tau2
	static sse(xs, ys, beta, tau1, tau2) {
		let s = 0;
		for (let i = 0; i < xs.length; i++) {
			const b = NSSFit.basis(xs[i], tau1, tau2);
			const yhat = beta[0]*b[0] + beta[1]*b[1] + beta[2]*b[2] + beta[3]*b[3];
			const e = ys[i] - yhat;
			s += e * e;
		}
		return s;
	}

	// Fit beta for fixed taus
	static fitBetaForTaus(YearsArray, YieldArray, tau1, tau2) {
		const { xs, ys } = NSSFit.buildSamples(YearsArray, YieldArray);
		if (xs.length < 4) return { ok: false, reason: "Need at least 4 valid points.", xs, ys };

		const A = xs.map(t => NSSFit.basis(t, tau1, tau2));
		const sol = NSSFit.solveBetaLeastSquares(A, ys);
		if (!sol.ok) return { ok: false, reason: "Linear solve failed (singular).", xs, ys };

		const beta = sol.beta;
		const SSE = NSSFit.sse(xs, ys, beta, tau1, tau2);
		return { ok: true, beta, tau1, tau2, SSE, n: xs.length };
	}

	// Grid search taus, solve beta each time; return best
	static fit(YearsArray, YieldArray, opts = {}) {
		const tau1Min = opts.tau1Min ?? 0.05;
		const tau1Max = opts.tau1Max ?? 5.0;
		const tau1Step = opts.tau1Step ?? 0.05;

		const tau2Min = opts.tau2Min ?? 0.3;
		const tau2Max = opts.tau2Max ?? 15.0;
		const tau2Step = opts.tau2Step ?? 0.1;

		let best = null;

		for (let tau1 = tau1Min; tau1 <= tau1Max + 1e-12; tau1 += tau1Step) {
			for (let tau2 = tau2Min; tau2 <= tau2Max + 1e-12; tau2 += tau2Step) {
				// often constrain tau2 > tau1 to avoid redundancy
				if (tau2 <= tau1) continue;

				const r = NSSFit.fitBetaForTaus(YearsArray, YieldArray, tau1, tau2);
				if (!r.ok) continue;

				if (!best || r.SSE < best.SSE) best = r;
			}
		}

		return best ? { ok: true, ...best } : { ok: false, reason: "No fit found." };
	}

	// Single-point NSS evaluation (like polynomialEval)
	static evaluation(params, Year) {
		const { beta, tau1, tau2 } = params;
		const t = Year;
		const b = NSSFit.basis(t, tau1, tau2); // [1, f1, f2, f3]
		return beta[0]*b[0] + beta[1]*b[1] + beta[2]*b[2] + beta[3]*b[3];
	}
}

// ==============================================================================
//  NSS Continuous Fitter (Excel-like workflow)
//  - Start: beta=[0.01,0.01,0.01,0.01], lambda1=1, lambda2=1
//  - Objective: SSE = sum((y - yhat)^2), skip null/undefined
//  - Optimizer: Nelder–Mead (derivative-free continuous optimization)
// ==============================================================================
class NSSExcelSolver {
  // （あなたの buildSamples / nssYield / sse はそのまま）
  // nssYield 内で NSSFit.basis を使う前提 :contentReference[oaicite:2]{index=2}

  // -------- Numerical gradient (finite difference) --------
	// ---- helpers ----
	static _isFinite(x) { return Number.isFinite(x); }

	// ----- NSS basis functions (given t in years, tau1,tau2) -----
	static _factor(t, tau) {
		const x = t / tau;
		if (x === 0) return 1;
		return (1 - Math.exp(-x)) / x;
	}

	// Return basis vector b = [1, f1, f2, f3] so that y = beta · b
	static basis(t, tau1, tau2) {
		const f1 = NSSExcelSolver._factor(t, tau1);
		const e1 = Math.exp(-t / tau1);
		const f2 = f1 - e1;

		const f1b = NSSExcelSolver._factor(t, tau2);
		const e2 = Math.exp(-t / tau2);
		const f3 = f1b - e2;

		return [1, f1, f2, f3];
	}

	// If yields look like "percent points" (e.g. 3.44, 6.0) convert to decimals (0.0344, 0.06)
	// If yields already decimals (<= 0.5 usually), keep.
	static normalizeYieldScale(yArray) {
		const vals = yArray.filter(v => v != null && Number.isFinite(v));
		if (vals.length === 0) return { scale: 1, ys: yArray.slice() };
		const maxAbs = Math.max(...vals.map(v => Math.abs(v)));

		// Heuristic:
		// - If max > 1.5, it's likely "percent points" like 3.44 (meaning 3.44%)
		// - Excel percent-formatted cells typically store decimals (0.0344). We want THAT.
		if (maxAbs > 1.5) {
			return { scale: 0.01, ys: yArray.map(v => (v == null ? null : v * 0.01)) };
		}
		return { scale: 1, ys: yArray.slice() };
	}

	static buildSamples(yearArray, yArray) {
		const xs = [];
		const ys = [];
		for (let i = 0; i < yearArray.length; i++) {
			const t = yearArray[i];
			const y = yArray[i];
			if (!Number.isFinite(t) || t <= 0) continue;
			if (y == null || !Number.isFinite(y)) continue;
			xs.push(t);
			ys.push(y);
		}
		return { xs, ys };
	}

	// Evaluate NSS using your basis: yhat = beta · basis(t, lambda1, lambda2)
	static nssYield(t, beta, lambda1, lambda2) {
		const b = NSSExcelSolver.basis(t, lambda1, lambda2); // <- your existing function
		return beta[0]*b[0] + beta[1]*b[1] + beta[2]*b[2] + beta[3]*b[3];
	}

	static sse(xs, ys, beta, lambda1, lambda2) {
		// enforce positive lambdas (Excel also typically constrains via initial values / solver behavior)
		if (!(lambda1 > 0) || !(lambda2 > 0)) return Number.POSITIVE_INFINITY;

		let sum = 0;
		for (let i = 0; i < xs.length; i++) {
			const yhat = this.nssYield(xs[i], beta, lambda1, lambda2);
			const e = ys[i] - yhat;
			sum += e * e;
		}
		return sum;
	}

	// ---------------- Nelder–Mead ----------------
	// Minimizes f(x). Deterministic. Good for matching Excel's local optimum from given start.
	static nelderMead(f, x0, opts = {}) {
		const n = x0.length;
		const maxIter = opts.maxIter ?? 5000;
		const tolX = opts.tolX ?? 1e-12;
		const tolF = opts.tolF ?? 1e-14;

		const alpha = opts.alpha ?? 1;     // reflection
		const gamma = opts.gamma ?? 2;     // expansion
		const rho   = opts.rho   ?? 0.5;   // contraction
		const sigma = opts.sigma ?? 0.5;   // shrink

		// initial simplex around x0
		const step = opts.step ?? (() => {
			const s = Array(n).fill(0);
			for (let i = 0; i < n; i++) s[i] = (Math.abs(x0[i]) + 1) * 0.05;
			return s;
		})();

		const simplex = [];
		simplex.push({ x: x0.slice(), fx: f(x0) });

		for (let i = 0; i < n; i++) {
			const xi = x0.slice();
			xi[i] += Array.isArray(step) ? step[i] : step;
			simplex.push({ x: xi, fx: f(xi) });
		}

		const sortSimplex = () => simplex.sort((a, b) => a.fx - b.fx);

		const centroid = (excludeLast = true) => {
			const m = excludeLast ? n : n + 1;
			const c = Array(n).fill(0);
			for (let i = 0; i < m; i++) {
				const x = simplex[i].x;
				for (let k = 0; k < n; k++) c[k] += x[k];
			}
			for (let k = 0; k < n; k++) c[k] /= m;
			return c;
		};

		const add = (a, b, s = 1) => a.map((v, i) => v + s*b[i]);
		const sub = (a, b) => a.map((v, i) => v - b[i]);
		const mul = (a, s) => a.map(v => v * s);

		sortSimplex();

		for (let iter = 0; iter < maxIter; iter++) {
			sortSimplex();
			const best = simplex[0];
			const worst = simplex[n];
			const secondWorst = simplex[n - 1];

			// termination checks
			// 1) function spread
			const fSpread = Math.abs(worst.fx - best.fx);
			// 2) simplex size
			let xSpread = 0;
			for (let i = 1; i <= n; i++) {
				for (let k = 0; k < n; k++) xSpread = Math.max(xSpread, Math.abs(simplex[i].x[k] - best.x[k]));
			}
			if (fSpread < tolF && xSpread < tolX) {
				return { ok: true, x: best.x, fx: best.fx, iter };
			}

			const c = centroid(true);

			// reflection: xr = c + alpha*(c - worst)
			const xr = add(c, sub(c, worst.x), alpha);
			const fr = f(xr);

			if (fr < best.fx) {
				// expansion: xe = c + gamma*(xr - c)
				const xe = add(c, sub(xr, c), gamma);
				const fe = f(xe);
				simplex[n] = (fe < fr) ? { x: xe, fx: fe } : { x: xr, fx: fr };
				continue;
			}

			if (fr < secondWorst.fx) {
				simplex[n] = { x: xr, fx: fr };
				continue;
			}

			// contraction
			let xc;
			if (fr < worst.fx) {
				// outside contraction: xc = c + rho*(xr - c)
				xc = add(c, sub(xr, c), rho);
			} else {
				// inside contraction: xc = c + rho*(worst - c)
				xc = add(c, sub(worst.x, c), rho);
			}
			const fc = f(xc);

			if (fc < worst.fx) {
				simplex[n] = { x: xc, fx: fc };
				continue;
			}

			// shrink towards best
			for (let i = 1; i <= n; i++) {
				const xs = add(best.x, sub(simplex[i].x, best.x), sigma);
				simplex[i] = { x: xs, fx: f(xs) };
			}
		}

		sortSimplex();
		return { ok: false, x: simplex[0].x, fx: simplex[0].fx, iter: maxIter };
	}

	// ---------------- Public API: fit like Excel ----------------
	// Input:
	// - yearArray: tenor in years (use full precision like Excel)
	// - yieldArray: observed yields (either decimals 0.0344 or percent points 3.44; auto-detect)
	static fit(yearArray, yieldArray, opts = {}) {
		const { ys: yNorm } = this.normalizeYieldScale(yieldArray);
		const { xs, ys } = this.buildSamples(yearArray, yNorm);

		if (xs.length < 4) {
			return { ok: false, reason: "Need at least 4 valid points (non-null yields)." };
		}

		// Excel start: betas 0.01, lambdas 1
		// To keep lambdas positive in an unconstrained optimizer, optimize u=ln(lambda).
		const x0 = [
			0.01, 0.01, 0.01, 0.01,
			Math.log(1), // u1
			Math.log(1), // u2
		];

		// objective on transformed vars
		const f = (x) => {
			const beta = [x[0], x[1], x[2], x[3]];
			const lambda1 = Math.exp(x[4]);
			const lambda2 = Math.exp(x[5]);
	
			return this.sse(xs, ys, beta, lambda1, lambda2);
		};

		const nm = this.bfgsMinimize(f, x0, {
			maxIter: opts.maxIter ?? 8000,
			tolX: opts.tolX ?? 1e-14,
			tolF: opts.tolF ?? 1e-16,
			step: opts.step ?? [0.02, 0.02, 0.02, 0.02, 0.5, 0.5], // important for converging to same local min
		});

		const xBest = nm.x;
		const beta = [xBest[0], xBest[1], xBest[2], xBest[3]];
		const lambda1 = Math.exp(xBest[4]);
		const lambda2 = Math.exp(xBest[5]);
		const SSE = this.sse(xs, ys, beta, lambda1, lambda2);

		return {
			ok: nm.ok,
			beta,
			lambda1,
			lambda2,
			SSE,
			iter: nm.iter,
			n: xs.length,
		};
	}

	static evaluation(params, Year) {
		const { beta, lambda1, lambda2 } = params;
		if (!beta || beta.length !== 4) return NaN;
		if (!(lambda1 > 0) || !(lambda2 > 0)) return NaN;
		//return this.nssYield(Year, beta, lambda1, lambda2);

		// yhat is in "decimal" if input yields were percent-points (due to normalizeYieldScale)
		const yhat = this.nssYield(Year, beta, lambda1, lambda2);

		// Convert output to percent-points (0.008 -> 0.8, 0.0344 -> 3.44)
		return yhat * 100;
	}

  static gradFD(f, x, eps = 1e-6) {
    const g = Array(x.length).fill(0);
    const fx = f(x);
    for (let i = 0; i < x.length; i++) {
      const h = eps * (Math.abs(x[i]) + 1);
      const xp = x.slice();
      const xm = x.slice();
      xp[i] += h;
      xm[i] -= h;
      const fp = f(xp);
      const fm = f(xm);
      g[i] = (fp - fm) / (2 * h);
    }
    return { g, fx };
  }

  // -------- Basic linear algebra helpers --------
  static dot(a, b) { let s = 0; for (let i=0;i<a.length;i++) s += a[i]*b[i]; return s; }
  static add(a, b, s=1) { return a.map((v,i)=> v + s*b[i]); }
  static sub(a, b) { return a.map((v,i)=> v - b[i]); }
  static normInf(a) { return Math.max(...a.map(v => Math.abs(v))); }

  static matVec(M, v) {
    const n = v.length;
    const out = Array(n).fill(0);
    for (let i=0;i<n;i++){
      let s=0;
      for (let j=0;j<n;j++) s += M[i][j]*v[j];
      out[i]=s;
    }
    return out;
  }

  static eye(n) {
    const I = Array.from({length:n}, (_,i)=>Array.from({length:n},(_,j)=>(i===j?1:0)));
    return I;
  }

  // BFGS update: H <- (I - rho*s*y^T) H (I - rho*y*s^T) + rho*s*s^T
  static bfgsUpdate(H, s, y) {
    const n = s.length;
    const ys = this.dot(y, s);
    if (!Number.isFinite(ys) || Math.abs(ys) < 1e-16) return H; // skip update
    const rho = 1 / ys;

    // Build (I - rho*s*y^T) and (I - rho*y*s^T)
    const A = this.eye(n);
    const B = this.eye(n);
    for (let i=0;i<n;i++){
      for (let j=0;j<n;j++){
        A[i][j] -= rho * s[i] * y[j];
        B[i][j] -= rho * y[i] * s[j];
      }
    }

    // temp = A * H
    const temp = Array.from({length:n}, ()=>Array(n).fill(0));
    for (let i=0;i<n;i++){
      for (let j=0;j<n;j++){
        let sum = 0;
        for (let k=0;k<n;k++) sum += A[i][k]*H[k][j];
        temp[i][j]=sum;
      }
    }

    // Hnew = temp * B
    const Hnew = Array.from({length:n}, ()=>Array(n).fill(0));
    for (let i=0;i<n;i++){
      for (let j=0;j<n;j++){
        let sum = 0;
        for (let k=0;k<n;k++) sum += temp[i][k]*B[k][j];
        Hnew[i][j]=sum;
      }
    }

    // + rho*s*s^T
    for (let i=0;i<n;i++){
      for (let j=0;j<n;j++){
        Hnew[i][j] += rho * s[i]*s[j];
      }
    }
    return Hnew;
  }

  // -------- Project to simple bounds (Excelの制約っぽいもの) --------
  static projectBounds(x, bounds) {
    if (!bounds) return x;
    const out = x.slice();
    for (let i=0;i<out.length;i++){
      const b = bounds[i];
      if (!b) continue;
      if (b.min != null && out[i] < b.min) out[i] = b.min;
      if (b.max != null && out[i] > b.max) out[i] = b.max;
    }
    return out;
  }

  // -------- GRG-like optimizer: BFGS + Armijo backtracking --------
  static bfgsMinimize(f, x0, opts = {}) {
    const maxIter = opts.maxIter ?? 5000;
    const gradTol = opts.gradTol ?? 1e-8;
    const stepTol = opts.stepTol ?? 1e-12;
    const fdEps   = opts.fdEps   ?? 1e-6;

    const bounds = opts.bounds ?? null;

    let x = this.projectBounds(x0, bounds);
    let H = this.eye(x.length);

    let { g, fx } = this.gradFD(f, x, fdEps);

    for (let iter=0; iter<maxIter; iter++) {
      const gInf = this.normInf(g);
      if (gInf < gradTol) return { ok:true, x, fx, iter, gradInf:gInf };

      // p = -H g
      const Hg = this.matVec(H, g);
      let p = Hg.map(v => -v);

      // If not a descent direction, fall back to steepest descent
      if (this.dot(p, g) >= 0) p = g.map(v => -v);

      // Armijo line search
      let alpha = 1.0;
      const c1 = 1e-4;
      const dg = this.dot(p, g); // negative expected

      let xNext = null, fNext = null;
      for (let ls=0; ls<50; ls++) {
        const trial = this.projectBounds(this.add(x, p, alpha), bounds);
        const ft = f(trial);
        if (Number.isFinite(ft) && ft <= fx + c1 * alpha * dg) {
          xNext = trial; fNext = ft;
          break;
        }
        alpha *= 0.5;
      }

      if (!xNext) return { ok:false, x, fx, iter, reason:"Line search failed" };

      const stepInf = this.normInf(this.sub(xNext, x));
      if (stepInf < stepTol) return { ok:true, x:xNext, fx:fNext, iter, gradInf:gInf };

      // Next gradient
      const { g: g2 } = this.gradFD(f, xNext, fdEps);

      // BFGS update
      const s = this.sub(xNext, x);
      const y = this.sub(g2, g);
      H = this.bfgsUpdate(H, s, y);

      x = xNext; fx = fNext; g = g2;
    }

    return { ok:false, x, fx, iter:maxIter, reason:"Max iterations" };
  }

  // ---------------- fit (GRG-like) ----------------
  static fitGRG(yearArray, yieldArray, opts = {}) {
    // 元の挙動（スケール推定）を維持するなら normalizeYieldScale をそのまま使う :contentReference[oaicite:3]{index=3}
    const { ys: yNorm } = this.normalizeYieldScale(yieldArray);
    const { xs, ys } = this.buildSamples(yearArray, yNorm);

    if (xs.length < 4) return { ok:false, reason:"Need at least 4 valid points." };

    // Excel start: beta=0.01, lambdas=1（あなたの最新仕様）
    const x0 = [0.01,0.01,0.01,0.01, Math.log(1), Math.log(1)];

    // Optional: Excelっぽく λ を現実的レンジに縛る（局所解暴走防止）
    // u=ln(lambda) なので bounds は log 空間で。
    const bounds = opts.bounds ?? [
      null, null, null, null,
      { min: Math.log(0.05), max: Math.log(50) }, // lambda1 in [0.05, 50]
      { min: Math.log(0.05), max: Math.log(50) }, // lambda2 in [0.05, 50]
    ];

    const f = (x) => {
      const beta = [x[0],x[1],x[2],x[3]];
      const lambda1 = Math.exp(x[4]);
      const lambda2 = Math.exp(x[5]);
      return this.sse(xs, ys, beta, lambda1, lambda2);
    };

    const sol = this.bfgsMinimize(f, x0, { ...opts, bounds });

    const xBest = sol.x;
    const beta = [xBest[0],xBest[1],xBest[2],xBest[3]];
    const lambda1 = Math.exp(xBest[4]);
    const lambda2 = Math.exp(xBest[5]);
    const SSE = this.sse(xs, ys, beta, lambda1, lambda2);

    return { ok: sol.ok, beta, lambda1, lambda2, SSE, iter: sol.iter, gradInf: sol.gradInf };
  }
}

class NSSExcelSolver2 {
	// ---- helpers ----
	static _isFinite(x) { return Number.isFinite(x); }

	// ----- NSS basis functions (given t in years, tau1,tau2) -----
	static _factor(t, tau) {
		const x = t / tau;
		if (x === 0) return 1;
		return (1 - Math.exp(-x)) / x;
	}

	// Return basis vector b = [1, f1, f2, f3] so that y = beta · b
	static basis(t, tau1, tau2) {
		const f1 = NSSExcelSolver._factor(t, tau1);
		const e1 = Math.exp(-t / tau1);
		const f2 = f1 - e1;

		const f1b = NSSExcelSolver._factor(t, tau2);
		const e2 = Math.exp(-t / tau2);
		const f3 = f1b - e2;

		return [1, f1, f2, f3];
	}

	// If yields look like "percent points" (e.g. 3.44, 6.0) convert to decimals (0.0344, 0.06)
	// If yields already decimals (<= 0.5 usually), keep.
	static normalizeYieldScale(yArray) {
		const vals = yArray.filter(v => v != null && Number.isFinite(v));
		if (vals.length === 0) return { scale: 1, ys: yArray.slice() };
		const maxAbs = Math.max(...vals.map(v => Math.abs(v)));

		// Heuristic:
		// - If max > 1.5, it's likely "percent points" like 3.44 (meaning 3.44%)
		// - Excel percent-formatted cells typically store decimals (0.0344). We want THAT.
		if (maxAbs > 1.5) {
			return { scale: 0.01, ys: yArray.map(v => (v == null ? null : v * 0.01)) };
		}
		return { scale: 1, ys: yArray.slice() };
	}

	static buildSamples(yearArray, yArray) {
		const xs = [];
		const ys = [];
		for (let i = 0; i < yearArray.length; i++) {
			const t = yearArray[i];
			const y = yArray[i];
			if (!Number.isFinite(t) || t <= 0) continue;
			if (y == null || !Number.isFinite(y)) continue;
			xs.push(t);
			ys.push(y);
		}
		return { xs, ys };
	}

	// Evaluate NSS using your basis: yhat = beta · basis(t, lambda1, lambda2)
	static nssYield(t, beta, lambda1, lambda2) {
		const b = NSSExcelSolver.basis(t, lambda1, lambda2); // <- your existing function
		return beta[0]*b[0] + beta[1]*b[1] + beta[2]*b[2] + beta[3]*b[3];
	}

	static sse(xs, ys, beta, lambda1, lambda2) {
		// enforce positive lambdas (Excel also typically constrains via initial values / solver behavior)
		if (!(lambda1 > 0) || !(lambda2 > 0)) return Number.POSITIVE_INFINITY;

		let sum = 0;
		for (let i = 0; i < xs.length; i++) {
			const yhat = this.nssYield(xs[i], beta, lambda1, lambda2);
			const e = ys[i] - yhat;
			sum += e * e;
		}
		return sum;
	}

	// ---------------- Nelder–Mead ----------------
	// Minimizes f(x). Deterministic. Good for matching Excel's local optimum from given start.
	static nelderMead(f, x0, opts = {}) {
		const n = x0.length;
		const maxIter = opts.maxIter ?? 5000;
		const tolX = opts.tolX ?? 1e-12;
		const tolF = opts.tolF ?? 1e-14;

		const alpha = opts.alpha ?? 1;     // reflection
		const gamma = opts.gamma ?? 2;     // expansion
		const rho   = opts.rho   ?? 0.5;   // contraction
		const sigma = opts.sigma ?? 0.5;   // shrink

		// initial simplex around x0
		const step = opts.step ?? (() => {
			const s = Array(n).fill(0);
			for (let i = 0; i < n; i++) s[i] = (Math.abs(x0[i]) + 1) * 0.05;
			return s;
		})();

		const simplex = [];
		simplex.push({ x: x0.slice(), fx: f(x0) });

		for (let i = 0; i < n; i++) {
			const xi = x0.slice();
			xi[i] += Array.isArray(step) ? step[i] : step;
			simplex.push({ x: xi, fx: f(xi) });
		}

		const sortSimplex = () => simplex.sort((a, b) => a.fx - b.fx);

		const centroid = (excludeLast = true) => {
			const m = excludeLast ? n : n + 1;
			const c = Array(n).fill(0);
			for (let i = 0; i < m; i++) {
				const x = simplex[i].x;
				for (let k = 0; k < n; k++) c[k] += x[k];
			}
			for (let k = 0; k < n; k++) c[k] /= m;
			return c;
		};

		const add = (a, b, s = 1) => a.map((v, i) => v + s*b[i]);
		const sub = (a, b) => a.map((v, i) => v - b[i]);
		const mul = (a, s) => a.map(v => v * s);

		sortSimplex();

		for (let iter = 0; iter < maxIter; iter++) {
			sortSimplex();
			const best = simplex[0];
			const worst = simplex[n];
			const secondWorst = simplex[n - 1];

			// termination checks
			// 1) function spread
			const fSpread = Math.abs(worst.fx - best.fx);
			// 2) simplex size
			let xSpread = 0;
			for (let i = 1; i <= n; i++) {
				for (let k = 0; k < n; k++) xSpread = Math.max(xSpread, Math.abs(simplex[i].x[k] - best.x[k]));
			}
			if (fSpread < tolF && xSpread < tolX) {
				return { ok: true, x: best.x, fx: best.fx, iter };
			}

			const c = centroid(true);

			// reflection: xr = c + alpha*(c - worst)
			const xr = add(c, sub(c, worst.x), alpha);
			const fr = f(xr);

			if (fr < best.fx) {
				// expansion: xe = c + gamma*(xr - c)
				const xe = add(c, sub(xr, c), gamma);
				const fe = f(xe);
				simplex[n] = (fe < fr) ? { x: xe, fx: fe } : { x: xr, fx: fr };
				continue;
			}

			if (fr < secondWorst.fx) {
				simplex[n] = { x: xr, fx: fr };
				continue;
			}

			// contraction
			let xc;
			if (fr < worst.fx) {
				// outside contraction: xc = c + rho*(xr - c)
				xc = add(c, sub(xr, c), rho);
			} else {
				// inside contraction: xc = c + rho*(worst - c)
				xc = add(c, sub(worst.x, c), rho);
			}
			const fc = f(xc);

			if (fc < worst.fx) {
				simplex[n] = { x: xc, fx: fc };
				continue;
			}

			// shrink towards best
			for (let i = 1; i <= n; i++) {
				const xs = add(best.x, sub(simplex[i].x, best.x), sigma);
				simplex[i] = { x: xs, fx: f(xs) };
			}
		}

		sortSimplex();
		return { ok: false, x: simplex[0].x, fx: simplex[0].fx, iter: maxIter };
	}

	// ---------------- Public API: fit like Excel ----------------
	// Input:
	// - yearArray: tenor in years (use full precision like Excel)
	// - yieldArray: observed yields (either decimals 0.0344 or percent points 3.44; auto-detect)
	static fit(yearArray, yieldArray, opts = {}) {
		const { ys: yNorm } = this.normalizeYieldScale(yieldArray);
		const { xs, ys } = this.buildSamples(yearArray, yNorm);

		if (xs.length < 4) {
			return { ok: false, reason: "Need at least 4 valid points (non-null yields)." };
		}

		// Excel start: betas 0.01, lambdas 1
		// To keep lambdas positive in an unconstrained optimizer, optimize u=ln(lambda).
		const x0 = [
			0.01, 0.01, 0.01, 0.01,
			Math.log(1), // u1
			Math.log(1), // u2
		];

		// objective on transformed vars
		const f = (x) => {
			const beta = [x[0], x[1], x[2], x[3]];
			const lambda1 = Math.exp(x[4]);
			const lambda2 = Math.exp(x[5]);
	
			return this.sse(xs, ys, beta, lambda1, lambda2);
		};

		const nm = this.nelderMead(f, x0, {
			maxIter: opts.maxIter ?? 8000,
			tolX: opts.tolX ?? 1e-14,
			tolF: opts.tolF ?? 1e-16,
			step: opts.step ?? [0.02, 0.02, 0.02, 0.02, 0.5, 0.5], // important for converging to same local min
		});

		const xBest = nm.x;
		const beta = [xBest[0], xBest[1], xBest[2], xBest[3]];
		const lambda1 = Math.exp(xBest[4]);
		const lambda2 = Math.exp(xBest[5]);
		const SSE = this.sse(xs, ys, beta, lambda1, lambda2);

		return {
			ok: nm.ok,
			beta,
			lambda1,
			lambda2,
			SSE,
			iter: nm.iter,
			n: xs.length,
		};
	}

	static evaluation(params, Year) {
		const { beta, lambda1, lambda2 } = params;
		if (!beta || beta.length !== 4) return NaN;
		if (!(lambda1 > 0) || !(lambda2 > 0)) return NaN;
		//return this.nssYield(Year, beta, lambda1, lambda2);

		// yhat is in "decimal" if input yields were percent-points (due to normalizeYieldScale)
		const yhat = this.nssYield(Year, beta, lambda1, lambda2);

		// Convert output to percent-points (0.008 -> 0.8, 0.0344 -> 3.44)
		return yhat * 100;
	}
}

class Regression {

	/**
	 * Polynomial regression (Least Squares Method)
	 * x: Array of x-values
	 * y: Array of y-values
	 * degree: Degree of the polynomial
	 * return: Coefficient array [a0, a1, a2, ...]
	 */
	static polynomialFit(x, y, degree) {
		// Clean null values
		const xs = [];
		const ys = [];

		for (let i = 0; i < x.length; i++) {
			if (y[i] != null && Number.isFinite(y[i])) {
				xs.push(x[i]);
				ys.push(y[i]);
			}
		}

		const n = degree + 1;
		const X = Array.from({ length: n }, () => Array(n).fill(0));
		const Y = Array(n).fill(0);

		for (let i = 0; i < xs.length; i++) {
			let xi = 1;
			const xiPowers = Array(n);
			for (let j = 0; j < n; j++) {
				xiPowers[j] = xi;
				xi *= xs[i];
			}

			for (let row = 0; row < n; row++) {
				for (let col = 0; col < n; col++) {
					X[row][col] += xiPowers[row] * xiPowers[col];
				}
				Y[row] += xiPowers[row] * ys[i];
			}
		}

		return this.gaussianElimination(X, Y);
	}

	/**
	 * Solves a system of linear equations using Gaussian elimination
	 * return: Solution vector
	 */
	static gaussianElimination(A, B) {
		const n = A.length;

		for (let i = 0; i < n; i++) {
			// Pivot selection
			let maxRow = i;
			for (let k = i + 1; k < n; k++) {
				if (Math.abs(A[k][i]) > Math.abs(A[maxRow][i])) {
					maxRow = k;
				}
			}

			// Swap rows
			[A[i], A[maxRow]] = [A[maxRow], A[i]];
			[B[i], B[maxRow]] = [B[maxRow], B[i]];

			// Normalize pivot row
			const div = A[i][i];
			for (let j = i; j < n; j++) A[i][j] /= div;
			B[i] /= div;

			// Eliminate other rows
			for (let k = 0; k < n; k++) {
				if (k === i) continue;
				const factor = A[k][i];
				for (let j = i; j < n; j++) {
					A[k][j] -= factor * A[i][j];
				}
				B[k] -= factor * B[i];
			}
		}

		return B;
	}

	/**
	 * Evaluates the polynomial using the given coefficients
	 * coef: Array [a0, a1, a2, ...]
	 * x: Input value
	 */
	static polynomialEval(coef, x) {
		return coef.reduce((sum, a, i) => sum + a * Math.pow(x, i), 0);
	}
}

// Fix dicimals at designated place
function fmtDecs(x, dc = 4, mode = "round") {
	if (!Number.isFinite(x)) return "—";

	let v = x;

	if (mode === "floor") {
		const f = 10 ** dc;
		v = Math.floor(x * f) / f;
	}
	else if (mode === "ceil") {
		const f = 10 ** dc;
		v = Math.ceil(x * f) / f;
	}
	// do nothing for "round"

	return v.toLocaleString('en-US', {
		minimumFractionDigits: dc,
		maximumFractionDigits: dc
	});
}

function NSS_equation(result) {
	if (!result.ok) return "Calculation Error";

	const beta = result.beta.map(b => b.toFixed(5));
	const lambda1 = result.lambda1.toFixed(5);
	const lambda2 = result.lambda2.toFixed(5);

	let lines = [];

	// β の表示
	lines.push(`β₀ = ${beta[0]}`);
	lines.push(`β₁ = ${beta[1]}`);
	lines.push(`β₂ = ${beta[2]}`);
	lines.push(`β₃ = ${beta[3]}`);

	// λ の表示
	lines.push(`λ₁ = ${lambda1}`);
	lines.push(`λ₂ = ${lambda2}`);

	return lines.join(", ");
}

function Poly_equation(coeffs) {
	const terms = [];
    coeffs.forEach((a, i) => {
        if (a === 0) return; // 0 の項は省略

        let sign = "";
        if (terms.length > 0) {
            sign = a >= 0 ? " + " : " - ";
        } else {
            sign = a >= 0 ? "" : "-";
        }

		const absA = Math.abs(a).toFixed(5);

        if (i === 0) {
            terms.push(`${sign}${absA}`);
        } else if (i === 1) {
            terms.push(`${sign}${absA}x`);
        } else {
            terms.push(`${sign}${absA}x^${i}`);
        }
	});
	return "y = " + terms.join("");
}
