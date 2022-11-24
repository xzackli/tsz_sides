using Pkg

Pkg.pkg"add WCS HDF5 Interpolations Cosmology JSON JLD2 PyPlot"
Pkg.pkg"add https://github.com/simonsobs/Pixell.jl.git"
Pkg.pkg"add https://github.com/xzackli/XGPaint.jl.git#prof"

using Downloads
Downloads.download("http://www.cita.utoronto.ca/~zack/cached_battaglia.jld2", 
    "cached_battaglia.jld2")
