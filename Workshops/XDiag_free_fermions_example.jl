using MKL
using LinearAlgebra,XDiag



println(BLAS.get_config())

#if you are getting some output like: 

#LBTConfig([ILP64] libmkl_rt.so, [LP64] libmkl_rt.so)

#then MKL is working


#after installing julia to try if the MKL and XDiag libraries work run the command:
# julia -t 16 --compile=all XDiag_free_fermions_example.jl
#here 16 is number of cores, put your own relevant to the machine used

# here the size is set so you have some time to see the utilisation of your cores,
# use task manager (windows) or htop (linux/mac) to see if all cores are utilised


let
    println("program starts")
    N = 26
    nfermions = 13

    
    block = Fermion(N, nfermions)
    ops = OpSum()
    for i in 1:N
        ops += "t" * Op("Hop", [i, mod1(i + 1, N)])
        #ops += "V" * Op("NN", [i, mod1(i + 1, N)])
    end
    ops["t"] = 1.0
    #ops["V"] = 2.0

    
    println("set hamiltonian, starting computation")
    e0 = eigval0(ops, block)
    println("ground state energy ",e0)
    println("computed")

end