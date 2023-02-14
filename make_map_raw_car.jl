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

box = [30   -30;           # RA
       -25     25] * Pixell.degree  # DEC
shape, wcs = geometry(CarClenshawCurtis{Float64}, box, 0.5 * Pixell.arcminute)

##
fid = h5open("/fs/lustre/cita/zack/projects/websky/big_box_m200c_v3.h5", "r")
ra, dec = deg2rad.(fid["ra"]), deg2rad.(fid["dec"])
redshift = collect(fid["redshift"])
halo_mass = collect(fid["halo_mass"])

##
perm = sortperm(dec, alg=ThreadsX.MergeSort)
ra = ra[perm]
dec = dec[perm]
redshift = redshift[perm]
halo_mass = halo_mass[perm]

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
function paint_map!(m, p::XGPaint.AbstractProfile, psa, sitp, masses, 
                    redshifts, αs, δs, irange)
    for i in irange
        α₀ = αs[i]
        δ₀ = δs[i]
        mh = masses[i]
        z = redshifts[i]
        θmax = XGPaint.θmax(p, mh * XGPaint.M_sun, z)
        profile_paint!(m, α₀, δ₀, psa, sitp, z, mh, θmax)
    end
end

function chunked_paint!(m, p::XGPaint.AbstractProfile, psa, sitp, masses, 
                        redshifts, αs, δs)
    m .= 0.0
    
    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize);
    
    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i
        i1, i2 = chunks[chunk_i]
        paint_map!(m, p, psa, sitp, masses, redshifts, αs, δs, i1:i2)
    end

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i - 1
        i1, i2 = chunks[chunk_i]
        paint_map!(m, p, psa, sitp, masses, redshifts, αs, δs, i1:i2)
    end
end

##
# cut = eachindex(halo_mass)[begin:end-5]
m = Enmap(zeros(shape), wcs)

print("Painting map.\n")
@time chunked_paint!(m, p, psa, sitp, halo_mass, redshift, ra, dec)


##
write_map(
    "map_raw.fits",
    m)

##
using PyPlot
plt.clf()
plt.figure()
plt.imshow(log10.(m.data'))
plt.axis("off")
plt.savefig("test.png", bbox_inches="tight",pad_inches = 0)
plt.gcf()



# ##
# m .= 0
# profile_paint!(m, 0., 0., psa, sitp, 0.1, 1e14, π/90)

# ##
# using Plots
# Plots.plot(log10.(m))
