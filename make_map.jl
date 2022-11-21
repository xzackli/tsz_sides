using Pixell, WCS, XGPaint
using Cosmology
using Interpolations
import XGPaint: AbstractProfile
using HDF5
import JSON
using JLD2

print("Threads: ", Threads.nthreads(), "\n")

stringdata=join(readlines("wcs.json"))
wcsdict=JSON.parse(stringdata)
shape = Int.((wcsdict["shape"][1], wcsdict["shape"][2]))
cdelt = Float64.((wcsdict["cdelt"][1], wcsdict["cdelt"][2]))
crpix = Float64.((wcsdict["crpix"][1], wcsdict["crpix"][2]))
crval = Float64.((wcsdict["crval"][1], wcsdict["crval"][2]))

wcs = Gnomonic{Float64}(
    cdelt,    # cdelt
    crpix,    # crpix
    crval,    # crval
    π/180     # "deg" to rad
)

##
fid = h5open("little_box.h5", "r")

ra, dec = deg2rad.(fid["ra"]), deg2rad.(fid["dec"])
redshift = collect(fid["redshift"])
halo_mass = collect(fid["halo_mass"])

# precomputed sky angles
α_map, δ_map = posmap(shape, wcs)
psa = (sin_α=sin.(α_map), cos_α=cos.(α_map), sin_δ=sin.(δ_map), cos_δ=cos.(δ_map))

##
print("Precomputing the model profile grid.\n")

# set up a profile to paint
p = XGPaint.BattagliaProfile(Omega_c=0.2589, Omega_b=0.0486, h=0.6774)

# beam stuff (not used in this particular script)
N_logθ = 512
rft = RadialFourierTransform(n=N_logθ, pad=256)

model_file::String = "cached_battaglia.jld2"
if isfile(model_file)
    print("Found cached Battaglia profile model. Loading from disk.\n")
    model = load(model_file)
    prof_logθs, prof_redshift, prof_logMs, prof_y = model["prof_logθs"], 
        model["prof_redshift"], model["prof_logMs"], model["prof_y"]
else
    print("Didn't find a cached profile model. Computing and saving.\n")
    logθ_min, logθ_max = log(minimum(rft.r)), log(maximum(rft.r))
    @time prof_logθs, prof_redshift, prof_logMs, prof_y = profile_grid(p; 
        N_logθ=N_logθ, logθ_min=logθ_min, logθ_max=logθ_max)
    save(model_file, Dict("prof_logθs"=>prof_logθs, 
        "prof_redshift"=>prof_redshift, "prof_logMs"=>prof_logMs, "prof_y"=>prof_y))
end


itp = Interpolations.interpolate(log.(prof_y), BSpline(Cubic(Line(OnGrid()))))
sitp = scale(itp, prof_logθs, prof_redshift, prof_logMs);

##
function paint_map!(m, p, psa, sitp, masses, redshifts, αs, δs)
    for i in axes(masses, 1)
        α₀ = αs[i]
        δ₀ = δs[i]
        Ms = masses[i]
        z = redshifts[i]
        profile_paint!(m, α₀, δ₀, p, psa, sitp, z, Ms)
    end
end

m = Enmap(zeros(shape), wcs)

print("Painting map.\n")
@time paint_map!(m, p, psa, sitp, halo_mass, redshift, ra, dec)

##
write_map(
    "map.fits",
    Enmap(
        m.data, 
        WCS.WCSTransform(2;
            cdelt = collect(cdelt),
            ctype = ["RA---TAN", "DEC--TAN"],
            crpix = collect(crpix),
            crval = collect(crval))))

##

using PyPlot
plt.clf()
plt.imshow((m), clim=(0.0, 1e-5))
plt.axis("off")
plt.tight_layout()
plt.savefig("test.png", bbox_inches="tight")
plt.gcf()
