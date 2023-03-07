# Installation instructions

We have not yet registered the module with the Julia package manager, meaning it has to be installed using a git repository URL. This can be done in the REPL as follows:

```julia
using Pkg

Pkg.add(PackageSpec(url="https://github.com/kopperud/Pesto.jl.git"))
```

Even though the module is not registered, the package manager (`Pkg`) will automatically resolve and install any necessary dependencies. Loading the module is done as follows:

```julia
using Pesto
```

Since Julia is a JIT (just-in-time) compiled language, any code must be compiled before it can be run, including modules. To save some time, there is also a pre-compiling step the first time a module is loaded. This means we have to wait a short while. Once the module is finished pre-compiling, you are now ready!
