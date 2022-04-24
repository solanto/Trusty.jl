# Trusty

Trusty aims to make basic deformation analyses of trusses subject to loads as convenient as possible. It even supports units via [Unitful](https://painterqubits.github.io/Unitful.jl/stable/)!

Defining a truss is as easy as:

```julia
Truss(1:2 .=> Section(29000000psi, 6inch^2) .=> 2:3,
      [0 1 2
       0 1 0]ft,
      ([1, 3], [1, 3]))
```

## Installation

```julia
] add Trusty
```

## Documentation

Currently, exported names are documented. You can access this documentation via your IDE or the `?` command in the REPL.

