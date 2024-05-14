"""
Abstract type representing a clock model.
"""
abstract type ClockModel end

"""
Mutable struct representing a strict clock model.

# Fields
- `clock_rate::Float64`: The rate of the clock.
"""
mutable struct StrictClock <: ClockModel
    clock_rate::Float64
end
