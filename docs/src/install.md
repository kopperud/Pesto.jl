# Installation

## Installing Julia

Julia is a high-level programming language similar to R, Matlab or Python. 
However, it is also a high-performance language.
Julia programs can be similarly fast to those written in compiled languages like [C/C++/Fortran](https://julialang.org/benchmarks/).
To install Julia, follow the instructions on the official [website](https://julialang.org/downloads/).

## Editors

There are several options for how to work with Julia:

* **(Recommended)** [Visual Studio Code](https://code.visualstudio.com) is an integrated developer environment (IDE), which has both a file editor and a console for entering commands. If you have used RStudio before, Visual Studio Code will be very similar.
* You can also run Julia in a [Jupyter](http://jupyter.org) notebook, as is often done with Python projects.
* Alternatively, you can edit script files with your editor of choice (for example notepad/vim), and either copy-paste lines of copy into the Julia console or use the `include("script.jl")` command.

## Installing Pesto.jl (stable)
To install the latest release of `Pesto.jl`, enter this in Julia:
```julia
import Pkg
Pkg.add("Pesto")
```
The package manager (`Pkg`) will automatically resolve and install any necessary dependencies.

## Installing Pesto.jl (dev)
If you want to use the developmental version of `Pesto.jl`, you can install it using a git repository URL. This can be done as follows:

```julia
import Pkg
Pkg.add(url="https://github.com/kopperud/Pesto.jl")
```

## Loading Pesto.jl
Pesto can be loaded like so:
```julia
using Pesto
```

Since Julia is a JIT (just-in-time) compiled language, any code must be compiled before it can be run. To save some time, there is also a pre-compiling step the first time a module is loaded. This means we have to wait a short while. Once the module is finished pre-compiling, you are now ready!