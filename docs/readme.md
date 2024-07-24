## Build documentation

Instructions for how to build the documentation locally.

1. Start a `julia` REPL in the "./docs" directory
2. Activate the environment with `] activate .`
3. Run `include("make.jl.")`

This will run the build script. Assuming that there are no errors, static HTML files will be written to `./docs/build/`. To see the docs, you need to run a web server. For example:

```bash
julia -e 'using LiveServer; serve(dir=\"/path/to/Pesto.jl/docs/build\")'
``` 

