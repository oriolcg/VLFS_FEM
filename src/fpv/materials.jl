
using DrWatson
@quickactivate "Offshore PV"

export give_material

# function give_material(name::String)
#     Materials_dict = Dict( "steel" => [932.26032, 2.15e11, 0.2192 ], #store material properties: density (kg/m^3) and young's modulus N/m^2 and in case of hollow material ration Ix_hollow/Ix_solid
#                         "hdpe" => [250, 12e9,1], #Wrong values. Taken from Regeout thesis.
#                         "neopren" => [116, 560e3,1],
#                         "cfrp" => [960.432, 40e9, 0.8704], #Carbon fiber reinforced polymer based on https://doi.org/10.1016/j.matpr.2020.07.636, lower end is used (conservative)
#                         "gfrp" => [937.8, 35e9, 0.9375],  #Glass fiber
#                         "eps" => [15, 6.5e6, 1],    #
#                         "pvc" => [896.432, 52e6, 0.8704]
#     )

#     ρ = Materials_dict[name][1]
#     E = Materials_dict[name][2]
#     γ = Materials_dict[name][3]
    
#     return ρ,E,γ
#     #return keys(Materials)
# end


function give_material(name::String)
    Materials_dict = Dict( 
        "hdpePhil" => [250, 12e9, 1.0], #Wrong values. Taken from Regeout thesis.
        "hdpe" => [952, 500e6, 1.0], #solid
        "hdpe_E2p0" => [952, 1000e6, 1.0], #solid 2*E
        "hdpe_E0p5" => [952, 250e6, 1.0], #solid 0.5*E
        "hdpe_F050" => [952, 500e6, 0.50], #solid
        "hdpe_F025" => [952, 500e6, 0.75], #solid
        "cfrp" => [960.432, 40e9, 1.0], #Carbon fiber reinforced polymer based on https://doi.org/10.1016/j.matpr.2020.07.636, lower end is used (conservative)
    )

    ρ = Materials_dict[name][1]
    E = Materials_dict[name][2]
    hollowRatio = Materials_dict[name][3]
    
    return ρ, E, hollowRatio
    #return keys(Materials)
end




