using DrWatson
using Gridap
using DelimitedFiles
using DataFrames:DataFrame
using DataFrames:Matrix
using Revise
using Plots

function refCoeff(k, a1, a2, a3, dx12, dx13)
    
    β = dx12*k
    γ = dx13*k

    Dk = (sin.(β)).^2 + (sin.(γ)).^2 + (sin.(γ-β)).^2  
    Dk = 2*Dk

    R1k = (sin.(β)).^2 + (sin.(γ)).^2

    R2k = sin.(γ) .* sin.(γ-β)

    R3k = -sin.(β) .* sin.(γ-β)

    Q1k = sin.(β) .* cos.(β) + sin.(γ) .* cos.(γ)

    Q2k = sin.(γ) .* cos.(γ-β) - 2*sin.(β)

    Q3k = sin.(β) .* cos.(γ-β) - 2*sin.(γ)

    Zi =  ( (a1 .* (R1k + im*Q1k)) + 
            (a2 .* (R2k + im*Q2k)) + 
            (a3 .* (R3k + im*Q3k)) ) ./ Dk
    
    Zr =  ( (a1 .* (R1k - im*Q1k)) + 
            (a2 .* (R2k - im*Q2k)) + 
            (a3 .* (R3k - im*Q3k)) ) ./ Dk
   
    return (Zi, Zr)
end

caseName = "cnst2_tau0p0_Tension2"
wrkDir = "data/sims_202301/mem_freq_damp_spec_cnst_tau0p0/mem_freq_damp_spec_"*caseName*"/"
H0 = 10

da = wload(wrkDir*"mem_data.jld2")

prbDa = da["prbDa"]

η₀ = da["η₀"]

ω = da["ω"]

k = da["k"]

a1 = prbDa[:,6]
a2 = prbDa[:,7]
a3 = prbDa[:,8]
dx12 = 1.0
dx13 = 2.3

ra1 = conj(a1./abs.(a1))

a1 = a1 .* ra1
a2 = a2 .* ra1
a3 = a3 .* ra1

Zi, Zr = refCoeff(k, a1, a2, a3, dx12, dx13);

Zt = prbDa[:,16]

R = abs.(Zr) ./ abs.(Zi);

C_a = abs.(Zi).^2 - abs.(Zr).^2 - abs.(Zt).^2 
C_a = ifelse.(C_a .< 0, 0, C_a)
C_a = sqrt.(C_a) ./ abs.(Zi)

khmin = 0
khmax = 9

ymax = 1.1

plt1 = plot(k*H0, R, linewidth=3, 
    xlabel = "kh",
    ylabel = "Reflection", legend=false,
    title = "Reflection Coeff",
    xlims = (khmin, khmax),
    ylims = (0, ymax)  )      


plt2 = plot(k*H0, abs.(Zt) ./ abs.(Zi), linewidth=3, 
    xlabel = "kh",
    ylabel = "Transmission", legend=false,
    title = "Transmission Coeff",
    xlims = (khmin, khmax),
    ylims = (0, ymax)  )      

plt3 = plot(ω, R, linewidth=3, 
    xlabel = "ω (rad/s)",
    ylabel = "Reflection", legend=false,
    title = "Reflection Coeff",
    ylims = (0, ymax) )     


plt4 = plot(ω, abs.(Zt) ./ abs.(Zi), linewidth=3, 
    xlabel = "ω (rad/s)",
    ylabel = "Transmission", legend=false,
    title = "Transmission Coeff",
    ylims = (0, ymax) )

pltAll = plot(plt1, plt2, plt3, plt4, layout=4, dpi=330,
    plot_title=caseName, plot_titlefontsize=10)
savefig(pltAll,wrkDir*"analysis/reflection.png")


plt1 = plot(k*H0, C_a, linewidth=3, 
    xlabel = "kh",
    ylabel = "Absorb", legend=false,
    title = "Absorb Coeff",
    xlims = (khmin, khmax),
    ylims = (0, ymax)  )      


plt2 = plot(ω, C_a, linewidth=3, 
    xlabel = "ω (rad/s)",
    ylabel = "Absorb", legend=false,
    title = "Absorb Coeff",
    ylims = (0, ymax) )     


plt3 = plot(k*H0, abs.(Zi), linewidth=3, 
    xlabel = "kh",
    ylabel = "Incident", legend=false,
    title = "Calc. Incid.",
    xlims = (khmin, khmax)  )      


plt4 = plot(ω, abs.(Zi), linewidth=3, 
    xlabel = "ω (rad/s)",
    ylabel = "Incident", legend=false,
    title = "Calc. Incident")     
    

pltAll = plot(plt1, plt2, plt3, plt4, layout=4, dpi=330,
    plot_title=caseName, plot_titlefontsize=10)
savefig(pltAll,wrkDir*"analysis/absorb.png")


plt1 = plot(ω, η₀, linewidth=3, 
    xlabel = "ω (rad/s)",
    ylabel = "Amp (m)", legend=false,
    title = "Amplitude")
    

plt2 = plot(k*H0, η₀, linewidth=3, 
    xlabel = "kh",
    ylabel = "Amp (m)", legend=false,
    title = "Amplitude",
    xlims = (khmin, khmax)  )

plt3 = plot(k*H0, ω, linewidth=3, 
    ylabel = "ω (rad/s)",
    xlabel = "kh", legend=false,
    title = "Dispersion",
    xlims = (khmin, khmax)  )

plt4 = plot(k*H0, 2*π./ω, linewidth=3, 
    ylabel = "T (s)",
    xlabel = "kh", legend=false,
    title = "Dispersion",
    xlims = (khmin, khmax)  )      

pltAll = plot(plt1, plt2, plt3, plt4, layout=4, dpi=330,
    plot_title=caseName, plot_titlefontsize=10)
savefig(pltAll,wrkDir*"analysis/input.png")

