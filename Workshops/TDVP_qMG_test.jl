using ITensors, ITensorMPS,BenchmarkTools,Plots,LinearAlgebra,Polynomials,Dates,Base.Threads,ITensorCorrelators,LsqFit,FFTW,JLD2,Distributions,MKL
ENV["JULIA_NUM_THREADS"] = 16
println("Number of threads: ", nthreads())


include("functions_tdvp.jl")

function save_list_to_file(list, filename)
    open(filename, "w") do file
        for item in list
            println(file, item)
        end
    end
end


# XXZ chain (OBC); the bond (mid_n, mid_n+1) has coupling Jmid
function xxz(n,Deltaz,Jmid)
    os = OpSum()

    mid_n::Int64=n/2

    for j = 1:n-1
        J = (j == mid_n) ? Jmid : 1.0
        os +=  Deltaz*J, "Sz", j, "Sz", j+1
        os +=  J * 0.5, "S+", j, "S-", j+1
        os +=  J * 0.5, "S-", j, "S+", j+1
    end

    return os
end







Jmid_list=[0.0]     # middle bond before the quench (0.0 = two decoupled halves)
n_list=[48]

Deltaz_list=[1.0]

dt=0.1
nsteps=30

for Deltaz in Deltaz_list
    for n in n_list
        s = siteinds("S=1/2", n)
        for Jmid in Jmid_list
            println(Jmid)

            state2 = Vector{String}(undef, n)
            for n in 1:n
                state2[n] = "Dn"
            end
            for n in 1:n
                if isodd(n)
                    state2[n] = "Up"
                end
            end
            psi_init=  MPS(s, state2)

            # ground state with modified middle bond
            Hdmrg = MPO(xxz(n,Deltaz,Jmid), s)
            e0_dmrg, psi0dmrg = dmrg(Hdmrg, psi_init; nsweeps=10, maxdim=[20,50,100,200], cutoff=1e-10)
            @show inner(psi0dmrg', Hdmrg, psi0dmrg) / inner(psi0dmrg, psi0dmrg)

            # quench: homogeneous chain
            Htdvp = MPO(xxz(n,Deltaz,1.0), s)

            psi0=psi0dmrg
            compute_and_plot_entanglement_renyi(psi0,n,0.0,1,Jmid,Deltaz)

            for tim in 1:nsteps

                psi0 = tdvp(
                Htdvp,
                -im*dt,
                psi0;
                time_step=-im*dt,
                maxdim=1000,
                cutoff=1e-10,
                normalize=true,
                reverse_step=false,
                outputlevel=1,
                )
                @save "psi0_timestep_$(tim).jld2"
                compute_and_plot_entanglement_renyi(psi0,n,round(tim*dt,digits=6),1,Jmid,Deltaz)
            end


        end
    end
end
