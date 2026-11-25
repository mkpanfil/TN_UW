function create_and_save_plot(x, y, title_str, xlabel_str, linecolor, marker, filename,ymin,ymax)
    plot_obj = plot(x, y, legend=false,xlabel=xlabel_str, ylabel=title_str, linecolor=linecolor, marker=marker,margin=2Plots.mm, ylims=(ymin, ymax))
    savefig(plot_obj, filename)
    return plot_obj
end

function create_and_save_plot_2(x, y, title_str, xlabel_str, linecolor, marker, filename,ymin,ymax)
    plot_obj = plot(x, y, legend=false,xlabel=xlabel_str, ylabel=title_str, linecolor=linecolor, marker=marker,margin=2Plots.mm, ylims=(ymin, ymax), xlims=(0.0,max(x)*1.05))
    savefig(plot_obj, filename)
    return plot_obj
end

function check_type(x)
    println("The type of x is: ", typeof(x))
end
function save_list_to_file(list, filename)
    open(filename, "a") do file
        for item in list
            println(file, item)
        end
    end
end

function generate_list(num_elements,mid)
    result = [mid+2]
    current = mid+2
    step = 2
    
    for i in 1:num_elements-1
        current += step * (-1)^i
        push!(result, current)
        step += 2
    end
    
    return result
end


function create_neel_with_flips(sites, num_flips,H,nsweeps,maxdim,cutoff,noise,mid,L)
    # Initialize the state with "Dn" and set "Up" at odd positions
    state2 = Vector{String}(undef, L)
    for n in 1:L
        state2[n] = "Dn"
    end
    for n in 1:L
        if isodd(n)
            state2[n] = "Up"
        end
    end
    # Flip the selected number of sites
    if num_flips > 0
        mid = div(L, 2)
        indices = generate_list(num_flips,mid)
        for i in indices
            state2[i] = "Up"
        end
    end
    println("indices ",indices)

    # Create the MPS
    psi_init = MPS(sites, state2)
    en0=0
    en0, psi0 = dmrg(H, psi_init; nsweeps, maxdim, cutoff,noise)
    @show flux(psi0)
    return en0,psi0
end


function create_neel(sites, num_flips,H,nsweeps,maxdim,cutoff,noise,mid,L)
    # Initialize the state with "Dn" and set "Up" at odd positions
    state2 = Vector{String}(undef, L)
    for n in 1:L
        state2[n] = "Dn"
    end
    for n in 1:L
        if isodd(n)
            state2[n] = "Up"
        end
    end
    # Create the MPS
    psi_init = MPS(sites, state2)
    en0=0
    en0, psi0 = dmrg(H, psi_init; nsweeps, maxdim, cutoff,noise)
    @show flux(psi0)
    return en0,psi0
end


function ent(psi, b)
    psi = orthogonalize!(psi, b)
    U, S, V = svd(psi[b], (linkinds(psi, b-1)..., siteinds(psi, b)...))
    SvN = 0.0
    for n in 1:ITensors.dim(S, 1)
        p = S[n, n]^2
        if p > 0
            SvN -= p * log(p)
        end
    end
    return SvN
end

function ent_renyi(psi, b, alpha_renyi)
    psi = orthogonalize!(psi, b)
    U, S, V = svd(psi[b], (linkinds(psi, b-1)..., siteinds(psi, b)...))

    if alpha_renyi == 0
        tol=1e-12
        rank = 0
        for n in 1:ITensors.dim(S, 1)
            if S[n,n]^2 > tol
                rank += 1
            end
        end

        return log(rank)
    end

    if alpha_renyi == 1
        # fallback to von Neumann
        SvN = 0.0
        for n in 1:ITensors.dim(S, 1)
            p = S[n, n]^2
            if p > 0
                SvN -= p * log(p)
            end
        end
        return SvN
    end

    if alpha_renyi == "inf"
        pmax = 0.0
        for n in 1:ITensors.dim(S, 1)
            p = S[n,n]^2
            pmax = max(pmax, p)
        end

        return -log(pmax)
    end

    s = 0.0
    for n in 1:ITensors.dim(S, 1)
        p = S[n, n]^2
        if p > 0
            s += p^alpha_renyi
        end
    end

    return log(s) / (1 - alpha_renyi)
end


function ent_renyi0(psi, b)
    tol=1e-12
    psi = orthogonalize!(psi, b)
    _, S, _ = svd(psi[b], (linkinds(psi, b-1)..., siteinds(psi, b)...))

    rank = 0
    for n in 1:ITensors.dim(S, 1)
        if S[n,n]^2 > tol
            rank += 1
        end
    end

    return log(rank)
end

function ent_renyi_inf(psi, b)
    psi = orthogonalize!(psi, b)
    _, S, _ = svd(psi[b], (linkinds(psi, b-1)..., siteinds(psi, b)...))

    pmax = 0.0
    for n in 1:ITensors.dim(S, 1)
        p = S[n,n]^2
        pmax = max(pmax, p)
    end

    return -log(pmax)
end







function entanglement_gap(psi, b)
    psi = orthogonalize!(psi, b)
    U, S, V = svd(psi[b], (linkinds(psi, b-1)..., siteinds(psi, b)...))
    singular_values = diag(S)
    probabilities = singular_values.^2
    # Sort probabilities in descending order
    sorted_probs = sort(probabilities, rev=true)
    gap = sorted_probs[1] - sorted_probs[2]

    
    return gap
end




function central_c_sat_metod(psi)
    mid::Int64 = length(psi)/2
    return (ent(psi, mid-1)-ent(psi, mid))
end

function compute_and_plot_entanglement(psi,t,DeltaL0,DeltaR0,bf)

    N = length(psi)
    mid::Int64=N/2

    entropies = Float64[]

    bonds = collect(1:N-1)
    
    for b in bonds
        push!(entropies, ent(psi, b))
    end
    textfilename = "Entanglement_entropy_ent_N_$(N)_DeltaL0_$(DeltaL0)_DeltaR0_$(DeltaR0)_OBC_stot_bonds__$(bf).csv"
    open(textfilename, "a") do file

        row = join(["$(t)"; map(string, entropies)], ",")
        println(file, row)
    end
end


function compute_and_plot_entanglement_renyi(psi,n,t,alpha_renyi,J2,Deltaz)

    N = length(psi)
    mid::Int64=N/2

    entropies = Float64[]

    bonds = collect(1:N-1)
    
    for b in bonds
        push!(entropies, ent_renyi(psi, b, alpha_renyi))
    end

    textfilename = "Entanglement_entropy_renyi_$(alpha_renyi)_n_$(n)_J2_$(J2)_Deltaz_$(Deltaz)_OBC_stot_.csv"
    open(textfilename, "a") do file

        row = join(["$(t)"; map(string, entropies)], ",")
        println(file, row)
    end
end


function compute_and_plot_entanglement_renyi0(psi,t,DeltaL0,DeltaR0)

    N = length(psi)
    mid::Int64=N/2

    entropies = Float64[]

    bonds = collect(1:N-1)
    
    for b in bonds
        push!(entropies, ent_renyi0(psi, b))
    end

    textfilename = "Entanglement_entropy_renyi0_N_$(N)_DeltaL0_$(DeltaL0)_DeltaR0_$(DeltaR0)_OBC_stot.csv"
    open(textfilename, "a") do file

        row = join(["$(t)"; map(string, entropies)], ",")
        println(file, row)
    end
end

function compute_and_plot_entanglement_renyi_inf_(psi,t,DeltaL0,DeltaR0)

    N = length(psi)
    mid::Int64=N/2

    entropies = Float64[]

    bonds = collect(1:N-1)
    
    for b in bonds
        push!(entropies, ent_renyi_inf(psi, b))
    end

    textfilename = "Entanglement_entropy_renyi_inf__N_$(N)_DeltaL0_$(DeltaL0)_DeltaR0_$(DeltaR0)_OBC_stot.csv"
    open(textfilename, "a") do file

        row = join(["$(t)"; map(string, entropies)], ",")
        println(file, row)
    end
end




#=
function mident(psi)
    L = length(psi)
    b::Int64=L/2
    psi = orthogonalize!(psi, b)
    U,S,V = svd(psi[b], (linkinds(psi, b-1)..., siteinds(psi, b)...))
    SvN = 0.0
    for n=1:ITensors.dim(S, 1)
        p = S[n,n]^2
        SvN -= p * log(p)
    end
    return SvN
end
=#

function calculate_chiral_operator_z(psi,i)
    c2xy_iip1::Float64= real(correlation_matrix(psi, "Sx", "Sy"; sites=i:i+1)[2])
    c2yx_iip1::Float64= real(correlation_matrix(psi, "Sy", "Sx"; sites=i:i+1)[2])
    return (c2xy_iip1-c2yx_iip1)
end

function calculate_chiral_operator_y(psi,i)
    c2zx_iip1::Float64= real(correlation_matrix(psi, "Sz", "Sx"; sites=i:i+1)[2])
    c2xz_iip1::Float64= real(correlation_matrix(psi, "Sx", "Sz"; sites=i:i+1)[2])
    return (c2zx_iip1-c2xz_iip1)
end


#corr_mat_zz=correlation_matrix(psi, "Sz", "Sz")
