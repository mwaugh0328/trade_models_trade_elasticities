### Eaton and Kortum, Exactly identified

#### Overview of code and description
---

- ``est_proc_ek_over.m``	main driver file, computes data moments, then finds the theta to best fit, computes standard errors.

- ``est_fun_over.m``	key function which simulates the model and constructs model implied moments via simmulation.

- `sim_trade_pattern_ek.m`	Generates trade flows and prices associated with the Eaton and Kortum (2002) model of trade.

- `sim_trade_pattern_ek_ponly.m`	Grabs only the prices, not the trade flows.

- `gen_moments.m`	generates the moments. TODO: make more transparent and faster.

- ` standard_error.m`	routine to construct bootstrap standard errors.

- `gen_fake_date.m`	given estimated parameters, constructs a fake data set prices and trade flows.

- `monte_carlo_proc.m` essentially estimates parameter values given the fake data set above.

- `gravregasym_logd.m` does the gravity regression with log distance.

- `obs_char.m`	computes observable characteristics of the price and trade data.

- `thetaest_est_exact.m`	this computes the DATA moments.

- `thetaest_est_mod_D.m`	computes the MODEL moments. This focus on the "Dni" moments.

- `trade_add_error.m` adds log normally distributed error to trade flows.
