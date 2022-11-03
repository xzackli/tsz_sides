# for_abhi

1. To install packages, run `julia install.jl`
2. Use `get_wcs.py` to configure a pixel map size and WCS information. This writes `wcs.json`.
3. Unzip some halos with `little_box.zip`.
4. Make maps with `julia -t 32 make_map.jl`. Here the number 32 refers to the thread count, but this is only used for the model precomputation at the moment.

<img src="test.png" 
alt="test"/></a>
