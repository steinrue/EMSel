import numpy as np
from scipy.stats import norm, binom, beta
from scipy.optimize import fsolve
from emsel.emsel_util import get_uq_a_exps, generate_states_new
from sympy import symbols, sympify, lambdify, diff
from numba import njit
from scipy.optimize import minimize as spoptmin
from tqdm import tqdm

class HMM:
    def __init__(self, num_approx, Ne, s_init=[0., 0.], s_bound=1, init_cond="theta", hidden_interp = "chebyshev", **kwargs):
        assert len(s_init)==2
        assert num_approx >= 0
        self.ZERO_ALLELE_TOLERANCE = 0
        self.init_cond = init_cond
        self.N = num_approx
        self.hidden_interp = hidden_interp
        self.Ne = Ne

        self.s1_init = self.s1 = s_init[0]
        self.s2_init = self.s2 = s_init[1]

        self.s_bound = s_bound
        self.gs, self.bounds = generate_states_new(self.N, self.Ne, hidden_interp)
        self.qs = 1-self.gs
        self.gs_product = np.multiply.outer(self.gs, self.gs)

        self.integral_bounds = self.bounds.copy()
        self.integral_bounds[0] = -np.inf
        self.integral_bounds[-1] = np.inf

        #self.a_ij = P(S_t+1 = j|S_t = i)
        self.a = self.calc_transition_probs_old([self.s1, self.s2])
        assert np.all(np.isclose(np.sum(self.a, axis=1), 1))
        self.a_init = self.a.copy()
        if self.init_cond == "uniform" or self.init_cond == "data_mean":
            if hidden_interp == "chebyshev":
                pre_istate = np.diff(self.bounds)
                self.init_state = pre_istate/np.sum(pre_istate)
            else:
                self.init_state = np.zeros_like(self.gs) + 1/self.N
        elif self.init_cond == "delta":
            self.init_state = np.zeros_like(self.gs)
            self.init_state[np.clip(np.argmin(np.abs(self.gs-kwargs["p"])), 1, self.N-2)] = 1.
        elif self.init_cond == "beta":
            self.init_state = np.zeros_like(self.gs)
            beta_param = kwargs["beta_coef"]
            beta_distrib = beta(beta_param, beta_param)
            beta_pdf = beta_distrib.pdf(self.gs[1:-1])
            self.init_state[1:-1] = beta_pdf/np.sum(beta_pdf)
        elif self.init_cond == "spikeandslab":
            self.init_state = np.ones_like(self.gs)
            self.init_state /= np.sum(self.init_state)
            self.init_state *= 1-kwargs["spike_frac"]
            self.init_state[np.clip(np.argmin(np.abs(self.gs - kwargs["spike_loc"])), 1, self.N - 2)] += kwargs["spike_frac"]
        else:
            raise TypeError("Invalid initial condition specification!")

        assert np.isclose(np.sum(self.init_state), 1.)
        self.init_init_state = np.copy(self.init_state)

    def get_update_func(self, update_type, update_args):
        if update_type == "full":
            return self.update_s_2D, []

        if update_type == "add":
            return self.update_s_add, []

        if update_type == "het":
            return self.update_s_het, []

        if update_type == "over":
            return self.update_s_over, []

        if update_type == "under":
            return self.update_s_under, []

        if update_type == "rec":
            return self.update_s_rec, []

        if update_type == "dom":
            return self.update_s_linear, [1., 1.]

        if update_type == "neutral":
            return self.update_s_neutral, []

        if update_type == "arc":
            self.multiplier_est = 0
            self.sym_s1, self.sym_s2, self.sym_lm = symbols("s1 s2 lm")
            try:
                self.constraint_expr = sympify(update_args["constraint_expr"])
            except:
                raise ValueError("Update type = arc but no constraint given!")
            self.s1_deriv = diff(self.constraint_expr, self.sym_s1)
            self.s2_deriv = diff(self.constraint_expr, self.sym_s2)
            return self.update_s_arc, [self.constraint_expr, self.s1_deriv, self.s2_deriv]

        if update_type == "linear":
            try:
                return self.update_s_linear, [update_args["alpha"], update_args["beta"]]
            except:
                raise ValueError("Update type = general linear but no alpha, beta given!")
        raise TypeError("Invalid update type!")


    def get_init_update_func(self, init_update_type):
        if init_update_type == "beta":
            return self.update_init_beta, self.beta_to_state_func, 2

        if init_update_type == "delta":
            return self.update_init_delta, self.delta_to_state_func, 1

        if init_update_type == "baumwelch":
            return self.update_init_baumwelch, self.baumwelch_to_state_func, self.N

        if init_update_type == "fixed":
            self.init_update_type = None
            return self.update_init_null, self.null_to_state_func, 1

        raise ValueError(f"Invalid init update type specified: {init_update_type}")

    def calc_transition_probs_old(self, s_vector):
        s1 = s_vector[0]
        s2 = s_vector[1]
        p_primes = np.clip(self.gs + self.gs * self.qs * (s2 * self.gs + s1 * (1 - 2 * self.gs)), 0, 1)
        sigmas = np.sqrt(self.gs * self.qs / (2*self.Ne))
        a_one = np.zeros((1, self.gs.shape[0]))
        a_one[:, 0] = 1.
        a_all = np.concatenate((a_one, np.diff(norm.cdf(np.expand_dims(self.integral_bounds, axis=-1),
                                            p_primes[1:-1], sigmas[1:-1]), axis=0).T, a_one[:, ::-1]), axis=0)
        return a_all

    def clip_and_renormalize(self, matrix, val):
        if self.hidden_interp == "chebyshev":
            uniform_matrix = np.diff(self.bounds)
            uniform_matrix /= np.sum(uniform_matrix)
            if uniform_matrix[0] < val:
                raise ValueError
            scale_factor = val/uniform_matrix[0]
            int_matrix = (1-scale_factor)*(matrix/np.sum(matrix)) + scale_factor*uniform_matrix
            int_matrix /= np.sum(int_matrix)
        else:
            int_matrix = np.clip(matrix, val, None)
            int_matrix /= np.sum(int_matrix)
        return int_matrix

    @staticmethod
    @njit(fastmath=True)
    def forward_one_numba(init_state, trans_matrix, a_t_to_bpmf, bpmf, sample_locs):
        N = init_state.shape[0]
        T = a_t_to_bpmf.shape[0]
        alphas_hat = np.zeros((T, N))
        cs = np.ones((T))
        cs[0] = 1. / np.sum(init_state * bpmf[a_t_to_bpmf[0], :])
        alphas_hat[0, :] = cs[0] * init_state * bpmf[a_t_to_bpmf[0], :]
        for t in np.arange(1, T):
            alphas_tilde = bpmf[a_t_to_bpmf[t], :] * np.dot(alphas_hat[t-1, :], trans_matrix)
            cs[t] = 1. / np.sum(alphas_tilde)
            alphas_hat[t, :] = cs[t] * alphas_tilde
        return alphas_hat.T, cs

    @staticmethod
    @njit(fastmath=True)
    def backward_one_numba(trans_matrix, a_t_to_bpmf, bpmf, cs):
        N = trans_matrix.shape[0]
        T = a_t_to_bpmf.shape[0]
        betas_hat = np.zeros((T, N))
        betas_hat[-1, :] = cs[-1]
        for t in range(2, T+1):
            const_fact = betas_hat[-t+1, :]*bpmf[a_t_to_bpmf[-t+1], :]
            betas_hat[-t, :] = cs[-t] * np.dot(trans_matrix, const_fact)
        return betas_hat.T

    @staticmethod
    @njit(fastmath=True)
    def xis_one_numba(alphas, betas, trans_matrix, bpmf_vals):
        N = alphas.shape[0]
        xis = np.zeros((N, N))
        bbetas = bpmf_vals[:, 1:]*betas[:, 1:]
        for i in np.arange(N):
            xis[i, :] +=  trans_matrix[i, :] * np.dot(bbetas, alphas[i, :-1])
        return xis

    def forward_backward_forward_one(self):
        self.alphas_hat, self.cs = self.forward_one_numba(self.init_state, self.a, self.a_t_to_bpmf_idx, self.bpmf_new, self.sample_locs)
        assert np.all(np.isclose(np.sum(self.alphas_hat, axis=0), 1.))
        self.betas_hat = self.backward_one_numba(self.a, self.a_t_to_bpmf_idx, self.bpmf_new, self.cs)
        if self.update_type not in ["add", "neutral"]:
            bpmf_vals = self.bpmf_new[self.a_t_to_bpmf_idx[:], :].T
            self.xis = self.xis_one_numba(self.alphas_hat, self.betas_hat, self.a, bpmf_vals)
            assert np.isclose(np.sum(self.xis), self.T-1)
        self.gammas = self.alphas_hat*self.betas_hat/self.cs
        assert np.all(np.isclose(np.sum(self.gammas, axis=0), 1.))

    def update_s_2D(self):
        post_ptp0 = np.sum(self.gammas[:, -1]*self.gs) - np.sum(self.gammas[:, 0]*self.gs)
        post_p3q = np.sum(np.dot(self.gammas[:, :-1].T, (self.gs**3*self.qs)))
        post_pq = np.sum(np.dot(self.gammas[:, :-1].T, (self.gs*self.qs)))
        post_p2q = np.sum(np.dot(self.gammas[:, :-1].T, (self.gs**2*self.qs)))
        post_pqm2p = np.sum(np.dot(self.gammas[:, :-1].T, (self.gs*self.qs*(1-2*self.gs))))
        post_p2qm2p = np.sum(np.dot(self.gammas[:, :-1].T, (self.gs**2*self.qs*(1-2*self.gs))))
        post_p2 = np.sum(np.dot(self.gammas[:, :-1].T, (self.gs**2)))
        post_pp1 = np.sum(self.xis*self.gs_product)
        post_pdeltap = post_pp1-post_p2
        s1_num_new = post_ptp0 * post_p3q - post_pdeltap * post_p2q
        s2_num_new = post_pqm2p * post_pdeltap - post_ptp0 * post_p2qm2p
        denom_new = post_pq * post_p3q - post_p2q**2
        return s1_num_new/denom_new, s2_num_new/denom_new

    def update_s_add(self):
        post_ptp0 = np.sum(self.gammas[:, -1]*self.gs) - np.sum(self.gammas[:, 0]*self.gs)
        post_pq = np.sum(np.dot(self.gammas[:, :-1].T, (self.gs*self.qs)))
        half_s = post_ptp0/post_pq
        return half_s, 2*half_s

    def update_s_het(self):
        ptp0 = np.sum(self.gammas[:, -1]*self.gs) - np.sum(self.gammas[:, 0]*self.gs)
        other = 2*(np.einsum("ij, i->", self.gammas[:, :-1], self.gs**2)-np.sum(self.xis*self.gs_product))
        denom = np.einsum("ij, i ->", self.gammas[:, :-1], (1-2*self.gs)**2*self.gs*self.qs)
        return (ptp0+other)/denom, 0.

    def update_s_under(self):
        s1, _ = self.update_s_het()
        if s1 >= 0:
            s1 = 0.
        return s1, 0.

    def update_s_over(self):
        s1, _ = self.update_s_het()
        if s1 <= 0:
            s1 = 0.
        return s1, 0.

    def update_s_neutral(self):
        return 0., 0.

    def update_s_rec(self):
        denom = np.sum(np.dot(self.gammas[:, :-1].T, (self.gs**3*self.qs)))
        num = np.sum(self.xis*self.gs_product)-np.sum(np.dot(self.gammas[:, :-1].T, (self.gs**2)))
        return 0., num/denom

    def update_s_linear(self, alpha, beta):
        denom_1 = beta*np.einsum("ij, i -> ", self.gammas[:, :-1], (1-2*self.gs)*self.gs*self.qs*(beta*(1-2*self.gs)+alpha*self.gs))
        denom_2 = alpha*np.einsum("ij, i ->", self.gammas[:, :-1], self.gs**2*self.qs*(beta*(1-2*self.gs)+alpha*self.gs))
        num_1 = beta*(np.sum(self.gammas[:, -1]*self.gs) - np.sum(self.gammas[:, 0]*self.gs))
        num_2 = (2*beta-alpha)*np.einsum("ij, i ->",self.gammas[:, :-1], self.gs**2)
        num_3 = (alpha-2*beta)*np.sum(self.xis*self.gs_product)
        num = num_1+num_2+num_3
        return beta*num/(denom_1+denom_2), alpha*num/(denom_1+denom_2)

    def beta_opt_func(self, beta_params, gammas):
        beta_density = beta.logpdf(self.gs[1:-1], beta_params[0], beta_params[1])
        betaneginf = np.isneginf(beta_density)
        if np.any(betaneginf):
            print("negative inf in beta fitting!!")
            beta_density[betaneginf] = -1000000
        return -np.sum(gammas[1:-1]*beta_density)

    def update_init_beta(self, gammas, min_init_val):
        fit_alpha_2, fit_beta_2 = spoptmin(self.beta_opt_func, np.array(self.init_params)+1, args=gammas, bounds=((0, None), (0, None)), method="Nelder-Mead")["x"]
        init_state = self.beta_to_state_func((fit_alpha_2, fit_beta_2), min_init_val)
        return init_state, (fit_alpha_2, fit_beta_2)

    def beta_to_state_func(self, fit_params, min_init_val):
        fit_alpha, fit_beta = fit_params
        ubounds = beta.cdf(self.bounds[1:], fit_alpha, fit_beta)
        lbounds = beta.cdf(self.bounds[:-1], fit_alpha, fit_beta)
        state_vals = ubounds-lbounds
        return self.clip_and_renormalize(state_vals, min_init_val)


    def update_init_delta(self, gammas, min_init_val):
        p_max = np.sum(self.gs*gammas)
        init_state = self.delta_to_state_func(p_max, min_init_val)

        return init_state, (p_max)


    def delta_to_state_func(self, p, min_init_val):
        state_vals = np.zeros_like(self.gs)
        state_vals[np.argmin(np.abs(self.gs-p))] = 1.

        return self.clip_and_renormalize(state_vals, min_init_val)


    def update_init_baumwelch(self, gammas, min_init_val):
        init_state = self.baumwelch_to_state_func(gammas, min_init_val)
        return init_state, init_state


    def baumwelch_to_state_func(self, gammas, min_init_val):
        return self.clip_and_renormalize(gammas, min_init_val)

    def update_init_null(self, gammas, min_init_val):
        init_init_state = self.null_to_state_func(0, min_init_val)
        return init_init_state, 0


    def null_to_state_func(self, param, min_init_val):
        return self.clip_and_renormalize(self.init_init_state, min_init_val)

    def get_lambdified_update_function(self, s1_coeff_1, s2_coeff_1, s1_coeff_2, s2_coeff_2, indep_coeff_1, indep_coeff_2, lambda_coeff_1, lambda_coeff_2, expr_3):
        symb_expr_1 = s1_coeff_1*self.sym_s1 + s2_coeff_1*self.sym_s2 + indep_coeff_1 - lambda_coeff_1*self.sym_lm
        symb_expr_2 = s1_coeff_2*self.sym_s1 + s2_coeff_2*self.sym_s2 + indep_coeff_2 - lambda_coeff_2*self.sym_lm
        symb_expr_3 = expr_3
        symb_exprs = [symb_expr_1, symb_expr_2, symb_expr_3]
        return lambdify(((self.sym_s1, self.sym_s2, self.sym_lm),),symb_exprs, "numpy")

    def update_s_arc(self, constraint_expr, s1_deriv, s2_deriv):
        s1_coeff_1 = -np.einsum("ij, i ->", self.gammas[:, :-1], (1 - 2 * self.gs) ** 2*self.gs*self.qs)
        s2_coeff_1 = -np.einsum("ij, i ->", self.gammas[:, :-1], (1 - 2 * self.gs)*self.gs ** 2*self.qs)
        s1_coeff_2 = s2_coeff_1
        s2_coeff_2 = -np.einsum("ij, i->", self.gammas[:, :-1], self.gs ** 3*self.qs)
        indep_coeff_1 = np.sum(self.gammas[:, -1] * self.gs) - np.sum(self.gammas[:, 0] * self.gs) + 2 * np.einsum("ij, i ->",
                self.gammas[:, :-1],self.gs ** 2) - 2 * np.sum(self.xis*self.gs_product)
        indep_coeff_2 = np.sum(self.xis*self.gs_product) - np.einsum("ij, i ->", self.gammas[:, :-1],self.gs ** 2)
        eval_exprs = self.get_lambdified_update_function(s1_coeff_1, s2_coeff_1, s1_coeff_2, s2_coeff_2, indep_coeff_1, indep_coeff_2, s1_deriv, s2_deriv, constraint_expr)
        soln = fsolve(eval_exprs, np.array([self.s1, self.s2, self.multiplier_est]))
        self.multiplier_est = soln[2]
        return soln[0], soln[1]

    def update_internals_from_datum(self, obs_counts, nts, sample_times):
        assert len(obs_counts) == len(nts)
        assert len(nts) == len(sample_times)
        self.T = int(sample_times[-1] + 1)
        self.sample_locs = sample_times

        self.obs_counts = np.zeros(self.T)
        self.obs_counts[self.sample_locs] = obs_counts
        self.nts = np.zeros(self.T, dtype=int)
        self.nts[self.sample_locs] = nts

        # self.b.pmf(count[t], nts[t], gs[i]) = P(a_t|n_t, f_ti)!
        # self.bpmf[t, i, n] = P(a_n_t|n_n_t, f_n_ti)

        # self.bpmf_new[k/n, i] = P(a_t = k|n_t = n, f_i)
        self.nts_uq = np.unique(self.nts)
        if 0 not in self.nts_uq:
            self.nts_uq = np.concatenate(([0], self.nts_uq))
        # emission probability mass function (bpmf)
        # self.bpmf_new = np.zeros((np.sum(self.nts_uq)+self.nts_uq.shape[0], self.N))
        self.bpmf_a = np.zeros(np.sum(self.nts_uq) + self.nts_uq.shape[0])
        self.bpmf_n = np.zeros_like(self.bpmf_a)

        self.bpmf_idx = np.cumsum(self.nts_uq + 1)
        for nt_i, nt in enumerate(self.nts_uq[1:], 1):
            self.bpmf_a[self.bpmf_idx[nt_i - 1]:self.bpmf_idx[nt_i]] = np.arange(nt + 1)
            self.bpmf_n[self.bpmf_idx[nt_i - 1]:self.bpmf_idx[nt_i]] = nt

        self.b = binom
        # self.bpmf = self.b.pmf(np.broadcast_to(self.obs_counts[..., None], self.obs_counts.shape+(self.N,)), np.broadcast_to(self.nts[..., None], self.nts.shape+(self.N,)), np.broadcast_to(self.gs, (self.nloc, self.T, self.N))).transpose([1,2,0])
        self.bpmf_new = self.b.pmf(np.broadcast_to(self.bpmf_a[..., None], self.bpmf_a.shape + (self.N,)),
                                   np.broadcast_to(self.bpmf_n[..., None], self.bpmf_n.shape + (self.N,)),
                                   np.broadcast_to(self.gs, (self.bpmf_a.shape[0], self.N)))

        self.a_t_to_bpmf_idx = np.zeros_like(self.nts)
        for t in np.nonzero(self.nts)[0]:
            self.a_t_to_bpmf_idx[t] = self.obs_counts[t] + self.bpmf_idx[
                np.where(self.nts_uq == self.nts[t])[0][0] - 1]

    def update_externals_from_datum(self, obs_counts, nts, sample_times):
        assert len(obs_counts) == len(nts)
        assert len(nts) == len(sample_times)
        T = int(sample_times[-1] + 1)
        sample_locs = sample_times

        full_obs_counts = np.zeros(T)
        full_obs_counts[sample_locs] = obs_counts
        full_nts = np.zeros(T, dtype=int)
        full_nts[sample_locs] = nts

        nts_uq = np.unique(full_nts)
        if 0 not in nts_uq:
            nts_uq = np.concatenate(([0], nts_uq))
        # emission probability mass function (bpmf)
        # .bpmf_new = np.zeros((np.sum(.nts_uq)+.nts_uq.shape[0], .N))
        bpmf_a = np.zeros(np.sum(nts_uq) + nts_uq.shape[0])
        bpmf_n = np.zeros_like(bpmf_a)

        bpmf_idx = np.cumsum(nts_uq + 1)
        for nt_i, nt in enumerate(nts_uq[1:], 1):
            bpmf_a[bpmf_idx[nt_i - 1]:bpmf_idx[nt_i]] = np.arange(nt + 1)
            bpmf_n[bpmf_idx[nt_i - 1]:bpmf_idx[nt_i]] = nt

        b = binom
        # .bpmf = .b.pmf(np.broadcast_to(.obs_counts[..., None], .obs_counts.shape+(.N,)), np.broadcast_to(.nts[..., None], .nts.shape+(.N,)), np.broadcast_to(.gs, (.nloc, .T, .N))).transpose([1,2,0])
        bpmf_new = b.pmf(np.broadcast_to(bpmf_a[..., None], bpmf_a.shape + (self.N,)),
                                   np.broadcast_to(bpmf_n[..., None], bpmf_n.shape + (self.N,)),
                                   np.broadcast_to(self.gs, (bpmf_a.shape[0], self.N)))

        a_t_to_bpmf_idx = np.zeros_like(full_nts)
        for t in np.nonzero(full_nts)[0]:
            a_t_to_bpmf_idx[t] = full_obs_counts[t] + bpmf_idx[
                np.where(nts_uq == full_nts[t])[0][0] - 1]

        return a_t_to_bpmf_idx, bpmf_new


    def compute_one_s(self, loc, tol, max_iter, min_init_val=1e-8, min_ic=5, save_history=False):
        s_hist = np.zeros((max_iter+1, 2))
        s_hist[0, :] = (self.s1, self.s2)
        ll_hist = np.zeros(max_iter+1)
        itercount = 0
        min_itercount = min_ic
        init_params_cache = np.zeros((min_itercount+1, self.init_update_size))
        first_five_stop = False
        invalid_s_flag = False
        self.init_params = np.zeros(self.init_update_size)
        while itercount < max_iter:
            self.forward_backward_forward_one()
            s1_next, s2_next = self.update_func(*self.update_func_args)
            ll_hist[itercount] = -np.sum(np.log(self.cs))
            if itercount > 0 and ll_hist[itercount] - ll_hist[itercount - 1] < tol:
                if itercount <= min_itercount:
                    first_five_stop = True
                elif first_five_stop:
                    best_idx = np.argmax(ll_hist[:min_itercount])
                    s_hist[best_idx+2:, :] = 0
                    ll_hist[best_idx+2:] = 0
                    ic_return = 12 if invalid_s_flag else 8
                    return s_hist if save_history else 0, [s_hist[best_idx+1, 0], s_hist[best_idx+1, 1]], ll_hist if save_history else ll_hist[ll_hist < 0][-1], init_params_cache[best_idx+1, :], best_idx, ic_return
                elif invalid_s_flag:
                    return s_hist if save_history else 0, [self.s1, self.s2], ll_hist if save_history else ll_hist[ll_hist < 0][-1], self.init_params, itercount, 4
                else:
                    return s_hist if save_history else 0, [self.s1, self.s2], ll_hist if save_history else ll_hist[ll_hist < 0][-1], self.init_params, itercount, 0
            if np.min(self.gammas[:, 0]) < -1e-3:
                print(f"idx {loc} iter {itercount} parameters absurd! {s1_next:.4f} {s2_next:.4f} {np.min(self.gammas[:, 0])}. Stopping.")
                return s_hist if save_history else 0, [self.s1, self.s2], ll_hist if save_history else ll_hist[ll_hist < 0][-1], self.init_params, itercount, 1
            if not self.is_valid_s(s1_next, s2_next):
                s_hist[itercount + 1, :] = self.s1, self.s2

                self.init_state, self.init_params = self.init_update_func(self.gammas[:, 0], min_init_val)
                if self.init_update_type is not None:
                    itercount += 1
                    invalid_s_flag = True
                    continue
                else:
                    return s_hist if save_history else 0, [self.s1, self.s2], ll_hist if save_history else ll_hist[ll_hist < 0][-1], self.init_params, itercount, 2
            if itercount > 1 and ll_hist[itercount] < ll_hist[itercount - 1] and ll_hist[itercount-1] < ll_hist[itercount-2]:
                if itercount <= min_itercount:
                    first_five_stop = True
                elif first_five_stop:
                    best_idx = np.argmax(ll_hist[:min_itercount])
                    s_hist[best_idx + 2:, :] = 0
                    ll_hist[best_idx + 2:] = 0
                    return s_hist if save_history else 0, s_hist[best_idx + 1, :], ll_hist if save_history else ll_hist[ll_hist < 0][-1], init_params_cache[best_idx+1, :], best_idx, 7
                else:
                    s_hist[itercount-1:, :] = 0
                    ll_hist[itercount-1:] = 0
                    return s_hist if save_history else 0, [s_hist[itercount-2, 0], s_hist[itercount-2, 1]], ll_hist if save_history else ll_hist[ll_hist < 0][-1], init_params_cache[-2, :], itercount-2, 3

            s_hist[itercount + 1, :] = (s1_next, s2_next)
            if self.update_type != "neutral":
                self.s1, self.s2 = s1_next, s2_next
                self.a = self.calc_transition_probs_old([self.s1, self.s2])

            init_params_cache[1:, :] = init_params_cache[:-1, :]
            init_params_cache[0, :] = self.init_params
            self.init_state, self.init_params = self.init_update_func(self.gammas[:, 0], min_init_val)
            assert np.all(np.isclose(np.sum(self.a, axis=1), 1))
            itercount += 1
        print(f"{loc} Max iterations exceeded without convergence. Stopping.")
        return s_hist if save_history else 0, [self.s1, self.s2], ll_hist if save_history else ll_hist[ll_hist < 0][-1], self.init_params, itercount-1, 9

    def compute_s(self, obs_counts, nts, sample_times, update_type, init_update_type, tol, max_iter, progressbar=True, data_mean=False, min_init_val=1e-8, save_history=False, **update_args):
        self.num_replicates = obs_counts.shape[0]
        self.update_type = update_type
        self.init_update_type = init_update_type
        self.update_func, self.update_func_args = self.get_update_func(update_type, update_args)
        self.init_update_func, self.init_params_to_state_func, self.init_update_size = self.get_init_update_func(init_update_type)
        self.s_history = np.zeros((max_iter+1, 2, self.num_replicates)) if save_history else 0
        self.ll_history = np.zeros((max_iter+1, self.num_replicates)) if save_history else 0
        self.s_final = np.zeros((2, self.num_replicates))
        self.ll_final = np.zeros((self.num_replicates, 1))
        self.init_params_array = np.zeros((self.num_replicates, self.init_update_size))
        self.itercount_array = np.zeros((self.num_replicates, 1), dtype=int)
        self.exit_codes = np.zeros((self.num_replicates))
        loop_idxs = tqdm(range(self.num_replicates)) if progressbar else range(self.num_replicates)
        for i in loop_idxs:
            self.update_internals_from_datum(obs_counts[i], nts[i], sample_times[i])
            if data_mean:
                mean_freq = np.sum(obs_counts[i]) / np.sum(nts[i])
                temp_init_state = np.zeros_like(self.init_state) + min_init_val
                temp_init_state[np.clip(np.argmin(np.abs(self.gs - mean_freq)), 1, self.N - 2)] = 1.
                self.init_state = temp_init_state / np.sum(temp_init_state)
            else:
                self.init_state = np.copy(self.init_init_state)
            self.s1 = self.s1_init
            self.s2 = self.s2_init
            self.a = self.a_init
            run_result = self.compute_one_s(i, tol, max_iter, min_init_val=min_init_val,save_history=save_history)
            self.s_final[:, i] = run_result[1]
            self.init_params_array[i, :] = run_result[3]
            self.itercount_array[i, 0] = run_result[4]
            self.exit_codes[i] = run_result[5]
            if save_history:
                self.s_history[:, :, i] = run_result[0]
                self.ll_history[:, i] = run_result[2]
                self.ll_final[i, 0] = run_result[2][run_result[2] < 0][-1]
            else:
                self.ll_final[i, 0] = run_result[2]
        return self.s_history, self.s_final, self.ll_history if save_history else self.ll_final, self.init_params_array, self.itercount_array, self.exit_codes

    def compute_multiple_ll(self, s1, s2, data_matrix, init_states=None):
        sample_locs_array = data_matrix[:, ::3]
        nts_array = data_matrix[:, 1::3]
        obs_counts_array = data_matrix[:, 2::3]
        ll_T = int(sample_locs_array[0, -1] + 1)
        ll_nloc = obs_counts_array.shape[0]
        ll_sample_locs = sample_locs_array[0, :]

        ll_obs_counts = np.zeros((ll_nloc, ll_T))
        ll_obs_counts[:, ll_sample_locs] = obs_counts_array
        ll_nts = np.zeros((ll_nloc, ll_T), dtype=int)
        ll_nts[:, ll_sample_locs] = nts_array

        # self.b.pmf(count[t], nts[t], gs[i]) = P(a_t|n_t, f_ti)!
        # self.bpmf[t, i, n] = P(a_n_t|n_n_t, f_n_ti)

        # self.bpmf_new[k/n, i] = P(a_t = k|n_t = n, f_i)
        ll_nts_uq = np.unique(ll_nts)

        # emission probability mass function (bpmf)
        # self.bpmf_new = np.zeros((np.sum(self.nts_uq)+self.nts_uq.shape[0], self.N))
        ll_bpmf_a = np.zeros(np.sum(ll_nts_uq) + ll_nts_uq.shape[0])
        ll_bpmf_n = np.zeros_like(ll_bpmf_a)

        ll_bpmf_idx = np.cumsum(ll_nts_uq + 1)
        for i, nt in enumerate(ll_nts_uq[1:], 1):
            ll_bpmf_a[ll_bpmf_idx[i - 1]: ll_bpmf_idx[i]] = np.arange(nt + 1)
            ll_bpmf_n[ll_bpmf_idx[i - 1]: ll_bpmf_idx[i]] = nt

        ll_b = binom
        # self.bpmf = self.b.pmf(np.broadcast_to(self.obs_counts[..., None], self.obs_counts.shape+(self.N,)), np.broadcast_to(self.nts[..., None], self.nts.shape+(self.N,)), np.broadcast_to(self.gs, (self.nloc, self.T, self.N))).transpose([1,2,0])
        ll_bpmf_new = ll_b.pmf(
            np.broadcast_to(ll_bpmf_a[..., None], (*ll_bpmf_a.shape, self.N)),
            np.broadcast_to(ll_bpmf_n[..., None], (*ll_bpmf_n.shape, self.N)),
            np.broadcast_to(self.gs, (ll_bpmf_a.shape[0], self.N)),
        )

        ll_a_t_to_bpmf_idx = np.zeros_like(ll_nts)
        for i, t in np.transpose(np.nonzero(ll_nts)):
            ll_a_t_to_bpmf_idx[i, t] = (
                    ll_obs_counts[i, t] + ll_bpmf_idx[np.where(ll_nts_uq == ll_nts[i, t])[0][0] - 1]
            )

        ll_a = self.calc_transition_probs_old([s1, s2])
        assert np.all(np.isclose(np.sum(ll_a, axis=1), 1))
        sample_times = np.nonzero(np.any(ll_nts, axis=0))[0]
        sample_time_diffs = np.diff(sample_times)
        uq_a_powers = np.unique(sample_time_diffs)
        uq_a_exps = get_uq_a_exps(ll_a, uq_a_powers)
        ll_cs = np.ones((ll_T, ll_nloc))
        if init_states is not None:
            assert init_states.shape == (ll_nloc, self.N)
            ll_alphas_tilde = np.einsum("ni, ni -> in", init_states, ll_bpmf_new[ll_a_t_to_bpmf_idx[:, 0], :])
        else:
            ll_alphas_tilde = np.einsum(
                "i, ni->in", self.init_state, ll_bpmf_new[ll_a_t_to_bpmf_idx[:, 0], :]
            )
        ll_cs[0, :] = 1.0 / np.sum(ll_alphas_tilde, axis=0)
        ll_alphas_hat = np.einsum("n, in -> in", ll_cs[0, :], ll_alphas_tilde)
        for i, t in enumerate(sample_times[1:]):
            ll_alphas_tilde = np.einsum(
                "in, ij, nj -> jn",
                ll_alphas_hat,
                uq_a_exps[np.where(uq_a_powers == sample_time_diffs[i])[0][0]],
                ll_bpmf_new[ll_a_t_to_bpmf_idx[:, t], :],
            )
            ll_cs[t, :] = 1.0 / np.sum(ll_alphas_tilde, axis=0)
            ll_alphas_hat = np.einsum("n, in -> in", ll_cs[t, :], ll_alphas_tilde)
            assert np.all(np.isclose(np.sum(ll_alphas_hat, axis=0), 1.0))
        return -np.sum(np.log(ll_cs), axis=0)

    def is_valid_s(self, s1, s2):
        if np.abs(s1) < self.s_bound and np.abs(s2) < self.s_bound:
            return True
        else:
            return False