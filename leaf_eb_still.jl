# This file declares all the necessary systems to recreate the energy balance model.

using Cropbox
using Unitful

"""
This system defines derived properties of the atmosphere for a given air temperature,
relative humidity, and pressure. Later systems use this as a mix-in to have access
to these parameters.
"""
@system Air(Controller) begin
    # Physical constants
    cp: heat_capacity_dry_air => 1010 ~ preserve(parameter, u"J / kg / K")
    Mair: molar_mass_dry_air => 0.029 ~ preserve(parameter, u"kg / mol")
    Mwater: molar_mass_water => 0.018 ~ preserve(parameter, u"kg / mol")
    Maw(Mair, Mwater): molar_mass_ratio => Mwater / Mair ~ track()
    R: gas_constant => 8.314 ~ preserve(parameter, u"J / K / mol")
    λ: heat_of_vaporization => 40.66e3 ~ preserve(parameter, u"J / mol")

    # Atmospheric conditions
    Ta: air_temperature => 298 ~ preserve(parameter, u"K")
    Pa: air_pressure => 100 ~ preserve(parameter, u"kPa")
    RH: relative_humidity => 0.5 ~ preserve(parameter)

    # Air densities
    rho(Pa, Mair, R, Ta): dry_air_density => begin
        Pa * Mair / (R * Ta) 
    end ~ track(u"kg / m^3")

    rho_mol(rho, Mair): dry_air_molar_density => begin
        rho / Mair
    end ~ track(u"mol / m^3")

    # Quantities related to water vapor
    γ(Pa, cp, Mwater, Maw, λ): psychrometric_constant => begin
        cp * Pa / (Maw * λ / Mwater)
    end ~ track(u"kPa / K")
    
    # https://en.wikipedia.org/wiki/Tetens_equation
    e_sat(Ta): saturation_vapor_pressure => begin
        k1 = 0.61078u"kPa"
        k2 = 17.27u"1"
        k3 = 237.3u"K"
        
        Ta_c = Ta - 273.15u"K"
        
        k1 * exp((k2 * Ta_c) / (Ta_c + k3))
    end ~ track(u"kPa")
    
    e(e_sat, RH): vapor_pressure => begin
        e_sat * RH
    end ~ track(u"kPa")

    vpd(e_sat, e): vapor_pressure_deficit => begin
        e_sat - e
    end ~ track(u"kPa")
    

    # Derivative of Tetens's formula
    d_esat(Ta): esat_slope => begin
        Ta_c = Ta - 273.15u"K"
        k1 = 0.61078u"kPa"
        k2 = 17.27u"1"
        k3 = 237.3u"K" 

        k1 * exp( (k2 * Ta_c) / (k3 + Ta_c) ) * k2 * k3 / (Ta_c + k3)^2
    end ~ track(u"kPa / K")

    ϵ(d_esat, γ): epsilon_ratio => begin
        d_esat / γ
    end ~ track() 
end

"""
This system defines some physical properties of the "leaf" that we are modeling.
These include absorption coefficients in the longwave and shortwave range
as well as the characteristic leaf dimension. The defaults correspond to a 
pine needle.
"""
@system Leaf(Controller) begin
    α_LW: emissivity => 0.98 ~ preserve(parameter)
    α_SW: absorptance => 0.50 ~ preserve(parameter)
    d: characteristic_dimension => 0.01 ~ preserve(parameter, u"m")
end

"""
The unified stomatal conductance model relates photosynthesis, vapor
pressure deficit, and atmospheric CO2 concentration with stomatal conductance.
Sensible defaults are chosen for a pine forest in 2015.
"""
@system UnifiedStomatalConductance(Controller, Air) begin
    g0: intercept => 0 ~ preserve(parameter, u"mmol / m^2 / s")
    # Parameter value is from an "arid pine site" in Metolius, OR
    g1: slope => 2.35 ~ preserve(parameter, u"kPa^0.5")

    # Paper used flux data from 2015 so Ca should match that time period
    Ca: co2_conc => 400 ~ preserve(parameter, u"μmol / mol")
    A: net_photosynthesis => 20 ~ preserve(parameter, u"μmol / m^2 / s")

    gs(g0, g1, A, Ca, vpd): stomatal_conductance => begin
       g0 + 1.6 * (1 + g1 / sqrt(vpd)) * A / Ca
    end ~ track(u"mmol / m^2 / s")
end

"""
This is a largely empirical model of conductances for a pine needle taken from the 
supplemental material of the source paper. This should be modified if studying
another system.
"""
@system NeedleleafConductance(Controller, Air, Leaf, UnifiedStomatalConductance) begin
    # Environmental input
    u: wind_speed => 0.5 ~ preserve(parameter, u"m/s")

    # Constants
    σ: stefan_boltzmann_const => 5.67e-8 ~ preserve(parameter, u"W / m^2 / K^4")

    # This is empirical so we have to try very hard to make Unitful happy
    gbH(u, d, rho_mol): boundary_layer_conductance => begin
        fudge_const = 0.006u"m^0.8 / s^0.4"
        fudge_const * rho_mol * u^0.6 / d^0.4
    end ~ track(u"mol / m^2 / s")

    gR(rho_mol, rho, α_LW, σ, Ta, cp): radiative_conductance => begin
        2 * rho_mol * (4 * α_LW * σ * Ta^3) / (rho * cp)
    end ~ track(u"mol / m^2 / s")

    gtot(gs, gbH): total_conductance => begin
        (gs^-1 + gbH^-1)^-1
    end ~ track(u"mol / m^2 / s")
end

"""
This system defines net radiation *at the leaf surface*. We can think of 
this as "isothermal" net radiation since it does not consider the temperature
of the leaf to calculate outgoing longwave.
"""
@system NetRadiation(Air, Leaf, Controller) begin
    # Environmental inputs
    SW_in:  shortwave_downwelling => 800 ~ preserve(parameter, u"W / m^2")
    SW_out: shortwave_upwelling => 100 ~ preserve(parameter, u"W / m^2")
    LW_in:  longwave_downwelling => 300 ~ preserve(parameter, u"W/ m^2")

    # Constants
    σ: stefan_boltzmann_const => 5.67e-8 ~ preserve(parameter, u"W / m^2 / K^4")

    SW_net(SW_in, SW_out, α_SW): shortwave_flux => begin
        (SW_in - SW_out) * α_SW
    end ~ track(u"W / m^2")

    LW_net(LW_in, Ta, α_LW, σ): longwave_flux => begin
        α_LW * (LW_in - σ * Ta^4)
    end ~ track(u"W / m^2")

    Rn(SW_net, LW_net): net_radiation => begin
        SW_net + LW_net
    end ~ track(u"W / m^2")
end

"""
The latent heat system models two endpoints: a "stomatally-imposed" transpiration
and an equilibrium, radiation-limited endpoint. Conductance values determine which
of these drive the overall latent heat flux through the decoupling coefficient Omega.
"""
@system LatentHeat(NeedleleafConductance, NetRadiation, Air) begin
    # First calculate the two endpoints for when omega = 0 or omega = 1
    LE_imp(gtot, λ, vpd, Pa): latent_heat_imposed => begin
        gtot * λ * vpd / Pa
    end ~ track(u"W / m^2")

    LE_eq(ϵ, Rn, gR, gbH): latent_heat_equilibrium => begin
        (ϵ * Rn) / (ϵ + 1 + gR / gbH)
    end ~ track(u"W / m^2")

    # Omega quantifies which endpoint is driving the observed flux
    Ω(ϵ, gR, gbH, gs): decoupling_coefficient => begin
        (ϵ + 2 + (gR / gbH)) / (ϵ + ((gR + gbH) / gs) + gR / gbH)
    end ~ track()

    LE(LE_imp, LE_eq, Ω): latent_heat => begin
        Ω * LE_eq + (1 - Ω) * LE_imp
    end ~ track(u"W / m^2")
end

"""
This final module actually calculates the delta-T term as specified in the paper. This system
should be manipulated when leaf temperature is the desired output.
"""
@system LeafTemperature_Still(LatentHeat, NeedleleafConductance, Air, Controller) begin
    Y(gbH, gR): conductance_ratio => gbH / (gbH + gR) ~ track()

    dT(Y, Rn, LE, cp, Mair, gbH): delta_t => begin
        Y * (Rn - LE) / (cp * Mair * gbH)
    end ~ track(u"K")

    T_leaf(dT, Ta): leaf_temperature => begin
        Ta + dT
    end ~ track(u"K")
end