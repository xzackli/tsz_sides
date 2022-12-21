using Pixell, WCS, XGPaint
using Cosmology
using Interpolations
import XGPaint: AbstractProfile
using HDF5
import JSON
using JLD2
using ThreadsX
using FileIO

print("Threads: ", Threads.nthreads(), "\n")

stringdata=join(readlines("wcs.json"))
wcsdict=JSON.parse(stringdata)

# In Python the shapes are backwards
shape = Int.((wcsdict["shape"][2], wcsdict["shape"][1]))
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


perm = sortperm(dec, alg=ThreadsX.MergeSort)
ra = ra[perm]
dec = dec[perm]
redshift = redshift[perm]
halo_mass = halo_mass[perm]

# precomputed sky angles
α_map, δ_map = posmap(shape, wcs)
psa = (sin_α=sin.(α_map), cos_α=cos.(α_map), sin_δ=sin.(δ_map), cos_δ=cos.(δ_map))

##
print("Getting the model profile grid.\n")

# set up a profile to paint
p = XGPaint.BattagliaProfile(Omega_c=0.2589, Omega_b=0.0486, h=0.6774)

# beam stuff (not used in this particular script)
rft = RadialFourierTransform(n=512, pad=256)


ells = rft.l
fwhm = deg2rad(8 * cdelt[1])  # heuristic: beam with FWHM of 8 × pixel size
σ² = fwhm^2 / (8log(2))
lbeam = @. exp(-ells * (ells+1) * σ² / 2)

##
model_file::String = "cached_battaglia_beamed.jld2"
if isfile(model_file)
    print("Found cached Battaglia profile model. Loading from disk.\n")
    model = load(model_file)
    logθs, redshifts, logMs, y_prof_grid = model["logθs"], 
        model["redshifts"], model["logMs"], model["y_prof_grid"]
else
    print("Didn't find a cached profile model. Computing and saving.\n")
    logθ_min, logθ_max = log(minimum(rft.r)), log(maximum(rft.r))
    @time logθs, redshifts, logMs, y_prof_grid = profile_grid(p; 
        N_logθ=length(rft.r), logθ_min=logθ_min, logθ_max=logθ_max);

    # now apply the beam
    @time XGPaint.transform_profile_grid!(y_prof_grid, rft, lbeam)
    y_prof_grid = y_prof_grid[begin+rft.pad:end-rft.pad, :, :]
    logθs = logθs[begin+rft.pad:end-rft.pad]
    XGPaint.cleanup_negatives!(y_prof_grid)

    jldsave(model_file; logθs, redshifts, logMs, y_prof_grid)
end

logθmax = XGPaint.build_max_paint_logradius(logθs, redshifts, logMs, y_prof_grid)
itp = Interpolations.interpolate(log.(y_prof_grid), BSpline(Cubic(Line(OnGrid()))))
sitp = scale(itp, logθs, redshifts, logMs);

##
function paint_map!(m, logθmax, psa, sitp, masses, 
                    redshifts, αs, δs, irange)
    for i in irange
        α₀ = αs[i]
        δ₀ = δs[i]
        mh = masses[i]
        z = redshifts[i]
        θmax = exp(logθmax(z, log10(mh)))
        profile_paint!(m, α₀, δ₀, psa, sitp, z, mh, θmax)
    end
end

function chunked_paint!(m, logθmax, psa, sitp, masses, 
                        redshifts, αs, δs)
    m .= 0.0
    
    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize);
    
    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i
        i1, i2 = chunks[chunk_i]
        paint_map!(m, logθmax, psa, sitp, masses, redshifts, αs, δs, i1:i2)
    end

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i - 1
        i1, i2 = chunks[chunk_i]
        paint_map!(m, logθmax, psa, sitp, masses, redshifts, αs, δs, i1:i2)
    end
end

m = Enmap(zeros(shape), wcs)

print("Painting map.\n")
@time chunked_paint!(m, logθmax, psa, sitp, halo_mass, redshift, ra, dec)

##
write_map(
    "map_beamed.fits",
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
plt.figure()
plt.imshow((m), clim=(0.0, 1e-5))
plt.axis("off")
plt.savefig("test_beamed.png", bbox_inches="tight",pad_inches = 0)
plt.gcf()
