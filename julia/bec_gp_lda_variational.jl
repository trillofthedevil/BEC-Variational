using Plots, Roots, Polynomials, LaTeXStrings, Integrals
using Optimization, OptimizationOptimJL

function psi(r,pvar)
  return pvar[1]*exp.(-pvar[2]*r.^2)
end

function dpsi(r,pvar)
  return -2*pvar[2]*r*psi(r,pvar)
end

function n0(r,pvar)
  return (psi(r,pvar)).^2 
end

function nex(n0as,N)
  return 8/(3*sqrt(pi))*(n0as)^1.5*sqrt(N)
end

function fNtot(pvar,model)
  A = pvar[1]; b = pvar[2]
  as = model[1]; N = model[2]
  sqpi = sqrt(pi)
  return pi*sqpi*A^2/b^1.5*(2^(-1.5)+8/3^2.5*A/sqpi*sqrt(N)*as^1.5)-1
end

function poly3(pvar,model)
  A = pvar[1]; b = pvar[2]
  as = model[1]; N = model[2]
  C = pi*sqrt(pi)/b^1.5
  return Polynomial( [-1, 0, C/2^1.5, C*8/3^2.5*sqrt(N/pi)*as^1.5] )
end

function f0(pvar,model)
  A = pvar[1]; b = pvar[2]
  as = model[1]; N = model[2]
  Ap = roots( poly3(pvar,model) )[end]
  n0 = 1/2^1.5
  return n0 / (n0+8/3^2.5*Ap*sqrt(N/pi)*as^1.5), Ap
end

function ntotr(r,pvar,model)
  n0r = n0(r,pvar)
  return n0r+nex.(n0r*model[1],model[2])
end 

function Etot(pvar,model)
  A = pvar[1]; b = pvar[2]
  as = model[1]; N = model[2]
  term1 = (3 * pi^1.5 * A^2) / (2^2.5 * b^0.5)
  term2 = (3 * pi^1.5 * A^2) / (2^4.5 * b^2.5)
  term3 = (pi^2.5 * as * N * A^4) / (4*b^1.5)
  term4 = (64 * pi^2 * N^1.5 * as^2.5 * A^5) / (3 * 5^1.5 * b^1.5)
  return (term1 + term2 + term3 + term4) / N
end

function cons(res, pvar, model) 
  A = pvar[1]; b = pvar[2]
  as = model[1]; N = model[2]
  res .= [fNtot(pvar,model)]
end

as = 5e-3
N = 1e6
Ap = 1.
b = 1.2

model = [as, N]
pvar = [Ap, b]

optprob = OptimizationFunction(Etot, Optimization.AutoForwardDiff(), cons = cons)
prob = OptimizationProblem(optprob, pvar, model, lcons = [0.], ucons = [0.])
sol = solve(prob, IPNewton())
pvar = sol.u

#sol = solve(prob, Ipopt.Optimizer())

println(sol, fNtot(pvar,model))

fc, Ap = f0(pvar,model)

N_int_sol = solve(IntegralProblem((r, p) -> 
                [4*pi*(r^2)*n0(r,p), 
                 4*pi*(r^2)*nex(n0(r,p)*model[1],model[2]), 
                 4*pi*(r^2)*ntotr(r,p,model)], 
                (0.0, Inf), pvar), QuadGKJL())

E_int_sol = solve(IntegralProblem((r, p) -> 
                [2*pi*(r^2)*dpsi(r,p)^2, 
                 2*pi*(r^4)*n0(r,p),
                 4*pi*(r^2)*2*pi*model[1]*model[2]*n0(r,p)^2,
                 4*pi*(r^2)*8*pi*model[1]*model[2]*n0(r,p)*nex(n0(r,p)*model[1],model[2])], 
                (0.0, Inf), pvar), QuadGKJL())

println(E_int_sol.u[1]," ", (3 * pi^1.5 * pvar[1]^2) / (2^2.5 * pvar[2]^0.5))
println(E_int_sol.u[2]," ", (3 * pi^1.5 * pvar[1]^2) / (2^4.5 * pvar[2]^2.5))
println(E_int_sol.u[3]," ", (pi^2.5 * model[1] * model[2] * pvar[1]^4) / (4*pvar[2]^1.5) )
println(E_int_sol.u[4]," ", (64 * pi^2 * model[2]^1.5 * model[1]^2.5 * pvar[1]^5) / (3 * 5^1.5 * pvar[2]^1.5) )

r = 0:0.01:15
n0r = n0(r,pvar) 
plot(r, n0r, line=(2,:dash), label=L"$n_c$")
plot!(r, nex.(n0r*model[1],model[2]), line=(2,:dot), label=L"$n_\mathrm{qd}$")
plot!(r, ntotr(r,pvar,model), line=(2,:solid), label=L"n_\mathrm{tot}")
plot!(xlabel=L"$r/a_0$", ylabel="density")
plot!(title=L"$f_0=%$(round(fc,digits=3)), \quad N_c, N_{qd}, N_{tot} = %$([round(u1,digits=3) for u1 in N_int_sol.u])$")
