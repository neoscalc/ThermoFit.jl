using Optim
using Plots

x0 = [0.5, 10]
norm = x0
x0 = x0 ./ norm
# define the rosenbrock function
function f(x, norm, a=0.1, b=0.02)
    x_denorm = x .* norm
    y = (a - x_denorm[1])^2 + b * (x_denorm[2] - x_denorm[1]^2)^2
    println("x = $x_denorm, y = $y")
    return y
end

# compute the function values on the grid
X = range(-1.0, 1.0, length=100)
Y = range(-1.0, 10.0, length=100)
X_norm = range(-1.0/norm[1], 1.0/norm[1], length=100)
Y_norm = range(-1.0/norm[2], 10.0/norm[2], length=100)
grid = Iterators.product(X_norm, Y_norm)
Z = [f(x, norm) for x in grid]

# plot contour lines of the function
contour(X, Y, Z', levels=150, c=:viridis)

# find mimimum of the function using the Nelder-Mead algorithm
res = optimize(x -> f(x, norm), x0, NelderMead())
scatter!([res.minimizer[1]]*norm[1], [res.minimizer[2]*norm[2]], label="minimum", c=:red)


