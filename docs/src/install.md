## Installation Julia

Go to the Julia website at [https://julialang.org/downloads/](https://julialang.org/downloads/), download the binary for your operating system, and save it somewhere appropriate. If you are using Mac/Linux, you may want to add the directory in which you saved the julia executable to the `PATH` variable. If you use the bash terminal, you can do this by appending the following to your `~/.bashrc` file:
```bash
export PATH="/path/to/julia-x.x.x/bin:$PATH"
```
Make sure to substitute `/path/to` with your directory, and fill in the correct version instead of `x.x.x`. This will allow you to start julia by typing `julia` into the terminal, regardless of which directory you are in.

## Installing Pesto.jl
We have not yet registered the module with the Julia package manager, meaning it has to be installed using a git repository URL. This can be done in the REPL as follows:

```julia
import Pkg

Pkg.add(url="https://github.com/kopperud/Pesto.jl")
```

Even though the module is not registered, the package manager (`Pkg`) will automatically resolve and install any necessary dependencies. Loading the module is done as follows:

```julia
using Pesto
```

Since Julia is a JIT (just-in-time) compiled language, any code must be compiled before it can be run. To save some time, there is also a pre-compiling step the first time a module is loaded. This means we have to wait a short while. Once the module is finished pre-compiling, you are now ready!
