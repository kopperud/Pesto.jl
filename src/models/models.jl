## abstract Model types

export Model

abstract type Model end

## Meaning only has one state
abstract type UniStateModel <: Model end 
## Meaning only has multiple states
abstract type MultiStateModel <: Model end


