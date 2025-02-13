#=
Here a series of pissoble loss functions should be defined.

All loss functions should be of the form:

function loss_function(x::Array{Float64,1}, y::Array{Float64,1})::Float64
    return loss::Float64
end

where x is the predicted value and y is the true value.


Couple of potentially useful losses...
for ASSEMBLAGE constraints
(1) Cross entropy

for COMPOSITION constraints
(2) q_cmp
(3) chi squared

=#


function chi_squared(y::AbstractVector, y_ref::AbstractVector)::Float64
    chi_squared = sum((y .- y_ref).^2 ./ y_ref)
    return chi_squared
end
