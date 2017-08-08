for f in ["basic","arpack1","medcplx","medrect","smrect","spdirect",
          "projection_and_transient","radius_abscissa","negproj","eigsnoconv",
          "numrange"
          ]
    println("Testing $f functionality")
    fnam = f*".jl"
    include(fnam)
end
