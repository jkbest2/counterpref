using AbstractGPs, Plots
using LinearAlgebra
using Random
using StatsBase

# Generate a GP realization
Random.seed!(2021)
n = 1000
X = range(0, 10, length = n)
# μf(x) = -0.1 * x * (x - 10)
μf(x) = 0.0
μ = μf.(X)
Σ = kernelmatrix(Matern32Kernel(), X)
U = cholesky(Σ)
Y = μ .+ U.L * randn(n)

plot(X, Y)

# samp = sample(1:length(X), weights(exp.(Y)), length(X); replace = false)
samp = reverse(sortperm(Y))
subsamp = samp[1:10:end]
# subsamp = reverse(sortperm(Y[1:10:end]))
# subsamp = reverse(sortperm(Y))
lowsamp = reverse(samp)
lowsub = lowsamp[1:10:end]

gp = GP(μf, Matern32Kernel())

plt = plot(X, Y; label = "True", color = "black")
savefig(plt, "figs/gp.svg")

nobs = 25
ridx = sample(1:length(X), nobs)
x = X[ridx]
y = Y[ridx]
fx = gp(x, 0.001)
p_fx = posterior(fx, y)
plt = plot(X, Y; label = "True", color = "black")
plot!(plt, p_fx, X; label = "Fit")
scatter!(plt, x, y; label = string(nobs) * " observations")
savefig(plt, "figs/" * string(nobs) * "_robs_gp.svg")

for nobs in [25, 50, 75, 100]
    x = X[subsamp[1:nobs]]
    y = Y[subsamp[1:nobs]]
    fx = gp(x, 0.001)
    p_fx = posterior(fx, y)
    plt = plot(X, Y; label = "True", color = "black")
    plot!(plt, p_fx, X; label = "Fit")
    scatter!(plt, x, y; label = string(nobs) * " observations")
    savefig(plt, "figs/" * string(nobs) * "_obs_gp.svg")
end

for ngrid in [5, 10, 15, 20]
    nobs = 25
    grid_idx = floor.(Int, range(1, n; length = ngrid))
    x = vcat(X[subsamp[1:nobs]],
            X[grid_idx])
    y = vcat(Y[subsamp[1:nobs]],
            Y[grid_idx])
    fx = gp(x, 0.001)
    p_fx = posterior(fx, y)
    plt = plot(X, Y; label = "True", color = "black")
    plot!(plt, p_fx, X; label = "Fit")
    scatter!(plt, x, y; label = string(nobs) * " observations")
    scatter!(plt, X[grid_idx], Y[grid_idx]; label = string(ngrid) * " grid obs")
    savefig(plt, "figs/" * string(ngrid) * "_grid_gp.svg")
end

function gap_idx(x; n = 1000)
    sort!(x)
    x = [1; x; n]
    d = diff(x)
    m, i = findmax(d)
    x[i] + m ÷ 2
end

function counter_idx(x, ncounter; n = 1000)
    cvec = [gap_idx(x)]
    for i in 2:ncounter
        append!(cvec, gap_idx([x; cvec]))
    end
    cvec
end

for ncounter in [5, 10, 15, 20]
    nobs = 25
    cidx = counter_idx(subsamp[1:nobs], ncounter)
    x = vcat(X[subsamp[1:nobs]],
            X[cidx])
    y = vcat(Y[subsamp[1:nobs]],
            Y[cidx])
    fx = gp(x, 0.001)
    p_fx = posterior(fx, y)
    plt = plot(X, Y; label = "True", color = "black")
    plot!(plt, p_fx, X; label = "Fit")
    scatter!(plt, x, y; label = string(nobs + ncounter) * " observations")
    scatter!(plt, X[cidx], Y[cidx]; label = string(ncounter) * " counter")
    savefig(plt, "figs/" * string(ncounter) * "_counter_gp.svg")
end
