using Dualization

import SCS
function scs_optimizer(;
    accel=10,
    alpha=1.5,
    eps=1e-9,
    max_iters=10_000,
    verbose=true,
)
    return JuMP.optimizer_with_attributes(
        SCS.Optimizer,
        "acceleration_lookback" => accel,
        "acceleration_interval" => 20,
        "alpha" => alpha,
        "eps_abs" => eps,
        "eps_rel" => eps,
        "max_iters" => max_iters,
        "verbose" => verbose,
    )
end
