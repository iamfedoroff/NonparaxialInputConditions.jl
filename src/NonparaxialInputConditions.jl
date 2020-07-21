module NonparaxialInputConditions

import ForwardDiff
import NLsolve

export nic_calculate


struct NIC{A<:AbstractArray}
    r :: A
    Y :: A
    D :: A
end


function nic_calculate(zs_func::Function, R::AbstractArray)

    function R_func(r)
        zs = zs_func(r)
        gamma = ForwardDiff.derivative(zs_func, r)
        return r + 2 * gamma / (1 - gamma^2) * zs
    end

    function D_func(r)
        zs = zs_func(r)
        gamma = ForwardDiff.derivative(zs_func, r)
        return 2 / (1 - gamma^2) * zs
    end

    function Y_func(r)
        R = R_func(r)
        dRdr = ForwardDiff.derivative(R_func, r)
        gamma = ForwardDiff.derivative(zs_func, r)
        return sqrt(r / R / dRdr * (1 + gamma^2) / (1 - gamma^2))
    end


    function func!(F, x)
        @. F = R_func(x) - R
        return nothing
    end

    sol = NLsolve.nlsolve(func!, R)


    r = sol.zero
    Y = @. Y_func(r)
    D = @. D_func(r)

    return NIC(r, Y, D)
end


end
