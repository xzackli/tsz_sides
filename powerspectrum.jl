using Pixell, Healpix
using FFTW
using PyPlot

extent_intermediate(shape, wcs) = (deg2rad.(shape .* wcs.cdelt))
area(shape, wcs) = abs(prod(extent_intermediate(shape, wcs)))
pixsize(shape, wcs) = area(shape, wcs) / prod(shape)
pixsize(m::Enmap{T, 2, A, Gnomonic{T}}) where {A,T} = m.wcs.cdelt
get_fft_norm(m::Enmap) = sqrt(pixsize(size(m), m.wcs) / prod(size(m)))

function laxes(shape, wcs)
    freqs = shape ./ extent_intermediate(shape, wcs)
    lx = fftfreq(shape[1], freqs[1]) .* 2π
    ly = fftfreq(shape[2], freqs[2]) .* 2π
    return lx, ly
end

function lmap(shape, wcs)
    lx, ly = laxes(shape, wcs)
    data = zeros(length(lx), length(ly), 2)
    for i in axes(data, 1)
        data[i,:,1] .= lx[i]
    end  
    for j in axes(data, 2)
        data[:,j,2] .= ly[j]
    end
    return Enmap(data, wcs)
end

function modlmap(shape, wcs)
    slmap = lmap(shape, wcs)
    l = sqrt.(dropdims(sum(slmap.^2; dims=3); dims=3))
    return Enmap(l, wcs)
end

function bin(bins, cl2d, wcs)
    lm = modlmap(size(cl2d), wcs)
    nbins = length(bins)
    counts = zeros(Int, nbins)
    ybins = zeros(nbins)
    for j in axes(cl2d,2), i in axes(cl2d, 1)
        ind = searchsortedfirst(bins, lm[i,j])
        if 0 < ind < nbins
            ybins[ind] += cl2d[i,j]
            counts[ind] += 1
        end
    end
    bin_centers = bins .+ (bins.step/2)
    return bin_centers, ybins ./ counts
end

##

m_websky = readMapFromFITS(joinpath("/fs/lustre/project/act/mocks/websky/v0.0/", 
    "tsz_8192.fits"), 1, Float64)
cl_websky = alm2cl(map2alm(m_websky))
ell_websky = 0:(length(cl_websky)-1)


##
m_raw = read_map("map_raw.fits"; trim=false)
# m_beamed = read_map("map_beamed.fits"; trim=false)

kmap_raw = fft(parent(m_raw)) .* get_fft_norm(m_raw)
# kmap_beamed = fft(parent(m_beamed)) .* get_fft_norm(m_beamed)

cl_2D_raw = real.(kmap_raw .* conj.(kmap_raw))
# cl_2D_beamed = real.(kmap_beamed .* conj.(kmap_beamed))

lmax = 10_000  
lb, cb_raw = bin(0:50:lmax, cl_2D_raw, m_raw.wcs)
# lb, cb_beamed = bin(0:50:lmax, cl_2D_beamed, m_beamed.wcs)
# the beam used in make_map_bandlimited.jl
# bl = gaussbeam(deg2rad(8 * m_raw.wcs.cdelt[1]), lmax + 1000);

##

plt.clf()
plt.plot(ell_websky, cl_websky .* ell_websky.^2 ./ 2π .* 1e12, lw=1, label="WebSky 0.4 nside 8192")
# plt.plot(lb, cb_beamed .* lb.^2 ./ bl[round.(Int,lb)].^2 ./ 2π .* 1e12, label="HD (realspace beam)")
plt.plot(lb, cb_raw .* lb.^2 ./ 2π .* 1e12, label="HD (raw)")
plt.xlim(100, lmax); plt.ylim(4e-2, 2)
plt.xscale("log"); plt.yscale("log")
plt.xlabel(raw"Multipole moment, $\ell$", weight="regular")
plt.ylabel(raw"$D_\ell$", weight="regular")
plt.legend()
plt.savefig("hirestest.png", dpi=200)
plt.gcf()

##

