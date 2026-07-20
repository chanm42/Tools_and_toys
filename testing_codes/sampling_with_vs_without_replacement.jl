using Statistics
using LinearAlgebra
using Random
using CairoMakie
using LaTeXStrings

# Software for simulating and plotting the difference between sampling with vs without replacement from a set
# The functions calc_sum_wr and calc_sum_wnr calculate the sum of a given set of integers 1:100 with and without replacement respectively
# This is a simple comparison of calculation of the averages of two vectors, but it demonstrates the power of sampling without replacement in convergence speed at the expense of more computational complexity / memory requirements

##================================================================##
function plot_things(x,y1,y2,y3,outfile)
##================================================================##
    # Creates the figure and sets the size
    f = Figure(size = (800,600))

    # Creates the axes on the figure, labels title and axes, sets gridlines to false and ticks to true for publication quality plots
    # tickalign is a value between 0 and 1, 1 is completely in, 0 is completely out
    ax = Axis(f[1,1], xlabel = "Number of Samples", ylabel = L"\bar{T}_{\mathrm{w}} \text{ vs } \bar{T}_{\mathrm{w/o}}",
        titlesize = 24,
        xlabelsize = 24, ylabelsize = 24,
        xlabelpadding = 3, ylabelpadding = 10,
        xgridvisible = true, ygridvisible = true,
        xminorgridvisible = false, yminorgridvisible = false,
        xticksvisible = true, yticksvisible = true,
        xticklabelsize = 20, yticklabelsize = 20,
        xminorticksvisible = true, yminorticksvisible = true,
        xticksmirrored = true, yticksmirrored = true,
        xtickalign = 0.66, ytickalign = 0.66,
        xminortickalign = 1, yminortickalign = 1,
        xticksize = 15, yticksize = 15,
        xminorticksize = 5, yminorticksize = 5,
        )
    # Copy and add scatterlines! commands as below for additional data vectors to be plotted
    # The marker keyword can differentiate between different data - see: https://docs.makie.org/stable/reference/plots/scatter#markers
    scatterlines!(ax, x, y1, color = :red, marker = :rect, label = L"\bar{T}_{\mathrm{w}}")
    scatterlines!(ax, x, y2, color = :blue, marker = :utriangle, label = L"\bar{T}_{\mathrm{w/o}}")
    lines!(ax, x, y3, color = :black, linestyle = :dash, label = L"\bar{T}_{\mathrm{exact}}")

    # Takes in x, yval, yval(err) as 3 separate vectors. Copy and add or remove lines from this for your data
    #errorbars!(ax, x, y1, err2, whiskerwidth = 10, color=:red)

    # Creates and edits the legend object - see: https://docs.makie.org/stable/reference/blocks/legend
    # Patchsize is the size of the legend rectangle, padding is the space between legend elements and the edges of the rectangle, margin is the space between the rectangle and axis lines
    # halign can be :left/:right/:center valign can be :top/:bottom/:center
    # labelsize is the font size of legend elements, markersize is the size of the markers
    leg = Legend(f[1,1], ax, patchsize = (20,20), padding = (10,10,10,10), margin = (20,20,20,20),
        halign = :right, valign = :bottom, tellwidth = false,
        labelsize = 24, markersize = 24)
    save(outfile,f)

    return nothing
end

function get_vec_swr(n)
    vec = zeros(n)
    for i=1:n
        vec[i] = i
    end
    return vec
end

function calc_sum_wr(num_samples,vec)
    n = length(vec)
    sum = 0.0e0
    for i = 1:num_samples
        r = rand(1:n)
        sum += vec[r]
    end
    sum /= num_samples
    return sum
end

function calc_sum_wnr(num_samples,vec)
    n = length(vec)
    rvec = shuffle(deepcopy(vec))
    sum = 0.0e0
    for i = 1:num_samples
        if i <= n
            x = popfirst!(rvec)
            sum += x
        end
    end
    num_samples <= n ? sum /= num_samples : sum /= n
    return sum
end

function calc_sum(num_samples,vec)
    sum_wr  = calc_sum_wr(num_samples,vec)
    sum_wnr = calc_sum_wnr(num_samples,vec)
    return sum_wr,sum_wnr
end

function main()
    n = 100
    vec  = get_vec_swr(n)
    sam_vec = []
    wr_vec  = []
    wnr_vec = []
    for num_samples = 10:10:1000
        sum_wr, sum_wnr = calc_sum(num_samples,vec)
        push!(sam_vec,num_samples)
        push!(wr_vec,sum_wr)
        push!(wnr_vec,sum_wnr)
    end
    exact_ans = sum(vec) / n
    exact_vec = fill(exact_ans,length(sam_vec))
    plot_things(sam_vec,wr_vec,wnr_vec,exact_vec,"swr_comp.pdf")
    return nothing
end

main()
