function placelegend()
    p = Plots.current()
    xl, yl = collect.(extrema.((Plots.xlims(p), Plots.ylims(p))))
    dx, dy = (xl[2] - xl[1])/3, (yl[2] - yl[1])/4
    x1, x2 = xl + [dx, -dx]
    y1, y2 = yl + [dy, -dy]
    tr = tl = br = bl = true
    for series in p.series_list
        x, y = series[:x], series[:y]
        tr && (tr = isnothing(findfirst(@. (x > x2) & (y > y2))))
        tl && (tl = isnothing(findfirst(@. (x < x1) & (y > y2))))
        br && (br = isnothing(findfirst(@. (x > x2) & (y < y1))))
        bl && (bl = isnothing(findfirst(@. (x < x1) & (y < y1))))
    end
    tr && return :topright
    tl && return :topleft
    br && return :bottomright
    bl && return :bottomleft
    return (0.5, 0.93)  # did not find empty corner, place legend outside
end