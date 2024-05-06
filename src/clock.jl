abstract type ClockModel end


mutable struct StrictClock <: ClockModel
    clock_rate::Float64
end