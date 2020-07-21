using NonparaxialInputConditions
using Test

f = 2.0

function zs_func(r)
    return r^2 / (4 * f)   # parabola
end

R = range(1e-3, 1.0, length=100)

rth = @. 2 * f^2 / R * (sqrt(1 + R^2 / f^2) - 1)
Yth = @. 1 - rth^2 / (4 * f^2)
Dth = @. 2 / (1 - rth^2 / (4 * f^2)) * rth^2 / (4 * f)

NIC = nic_calculate(zs_func, R)

@test isapprox(NIC.r, rth)
@test isapprox(NIC.Y, Yth)
@test isapprox(NIC.D, Dth)
