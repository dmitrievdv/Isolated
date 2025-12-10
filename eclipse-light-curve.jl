

# m_1 = 10
# m_2 = 15

# r = 0.1
# h = 0.5
# v = 0.3

# x_0 = √((r+1)^2 - h^2)

# xs = 1.5*collect(range(-x_0, x_0, 1000))
# ts = xs ./ v



function calc_intersect_area(h, x, r)
    s = sqrt(h^2 + x^2)
    if s ≤ 1 - r
        return 1.0
    end
    if s ≥ r + 1
        return 0.0
    end
    cosα = -(r^2 - 1 - s^2)/2/s
    cosβ = -(1^2 - r^2 - s^2)/2/s/r

    if abs(cosβ) > 1.0
        cosβ = cosβ/abs(cosβ)
    end
    # println("$h $x $r")
    α = acos(cosα)
    β = acos(cosβ)

    s_α = α - sin(α)*cos(α)
    s_β = β*r^2 - sin(β)*cos(β)*r^2

    return (s_α + s_β)/π/r^2
end

function calc_eclipse_magnitude(m_1, m_2, s)
    m_2_ecl = if 1 - s < 1e-6
        1000
    else
        m_2 - 2.5*log10(1 - s)
    end

    return -2.5*log10(10^(-0.4*m_1) + 10^(-0.4*m_2_ecl))
end

function calc_transit_mag(all_mag, primary_rel_flux, a, ϕ, ϕ₀, h, r)
    all_flux = 10^(0.4(10-all_mag))
    primary_flux = all_flux/(1 + 1/primary_rel_flux)
    secondary_flux = primary_flux/primary_rel_flux

    x = (ϕ - ϕ₀)*a
    s = calc_intersect_area(abs(h), x, abs(r))
    flux = all_flux - secondary_flux*s
    return 10 - 2.5*log10(flux)
end