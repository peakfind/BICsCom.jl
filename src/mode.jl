# functions for compute fields of modes corresponding to BICs

function compute_mode(sp::SamplingPoints, k, r, inn, ext, homo)
    n = sp.n
    cydc = build_cylinder_cache(4n, k, r, inn, ext)
    coef = get_coeff(4n, cydc, inn, ext)
    bydc = build_boundary_cache(sp, 4n, k, ext)
    Λ, A = assemble_dtn(sp, 4n, coef, bydc, k, ext)
    Δ = apply_bc(Λ, n)
    apply_tbc!(Δ, sp, k; homo = homo)
    
    evals, evecs = eigen(Δ)
    idx = argmin(abs.(evals))
    v = evecs[:, idx]
    
    @printf("k = %.4f, |λ_bic| = %.2e\n", k, abs(evals[idx]))
    
    # Reconstruct the Dirichlet data on 4 edges
    vb = v[1:n]
    vl = v[n+1:2n]
    vt = v[2n+1:end]
    vr = reverse(vl)
    u = [vb; vl; vt; vr]
    
    c = A \ u
    
    return c
end

function eval_field(xs, ys, k, r, inn, ext, dv, trunc, coeff)
     field = zeros(ComplexF64, length(ys), length(xs))
    ki = k * sqrt(inn)
    ke = k * sqrt(ext)
    
    orders = [j - 1 - trunc ÷ 2 for j in 1:trunc]

    for i in eachindex(xs), j in eachindex(ys)
        x = xs[i]
        y = ys[j]
        rd = sqrt(x^2 + y^2)
        θ = atan(y, x)
        
        val = zero(ComplexF64)
        
        if rd < r 
            for t in 1:trunc 
                val += dv[t] * besselj(orders[t], ki * rd) * exp(im * orders[t] * θ)
            end
        else
            for t in 1:trunc
                radial = coeff[1, t] * besselj(orders[t], ke * rd) +
                         coeff[2, t] * bessely(orders[t], ke * rd)
                val += dv[t] * radial * exp(im * orders[t] * θ)
            end
        end
        
        field[j, i] = val
        
    end
    
    return field
end