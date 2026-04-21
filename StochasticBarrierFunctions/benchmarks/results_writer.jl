using StochasticBarrierFunctions: eta, beta, psafe, total_time

const _CSV_HEADER = "tool,table,class,system,barrier_type,pwc_method,degree,num_partitions,time_s,eta,beta,Ps,status"

function results_dir()
    dir = get(ENV, "SB_RESULTS_DIR", joinpath(@__DIR__, "results"))
    mkpath(dir)
    return dir
end

function results_path()
    joinpath(results_dir(), "julia_results.csv")
end

function _ensure_header(path)
    if !isfile(path) || filesize(path) == 0
        open(path, "w") do io
            println(io, _CSV_HEADER)
        end
    end
end

# Parse yaml file path like
#   benchmarks/linear/systems/contraction/SOS/sos_deg2.yaml
#   benchmarks/linear/systems/contraction2/PWC/DUAL/contraction_dual_0.10.yaml
#   benchmarks/nonlinear/systems/pendulum/PWC/GD/pwc_gd_120.yaml
# Returns (table, class, system, barrier_type, pwc_method, degree_or_missing, num_partitions_or_missing)
function parse_yaml_metadata(yaml_file::AbstractString, config::AbstractDict)
    parts = splitpath(yaml_file)
    # locate "systems" segment
    idx = findfirst(==("systems"), parts)
    class_dir = idx === nothing ? "" : parts[idx - 1]
    system_dir = idx === nothing ? "" : parts[idx + 1]
    barrier_dir = idx === nothing ? "" : (length(parts) ≥ idx + 2 ? parts[idx + 2] : "")
    pwc_dir     = (barrier_dir == "PWC" && length(parts) ≥ idx + 3) ? parts[idx + 3] : ""

    system_map = Dict(
        "contraction"  => ("Linear",         "2D Contraction Map 1", "3"),
        "contraction2" => ("Linear",         "2D Contraction Map 2", "4"),
        "twotank"      => ("Linear",         "2D Two Tank",          "3"),
        "room"         => ("Linear",         "3D Room Temp.",        "3"),
        "quadrotor"    => ("Linear",         "6D Quadrotor",         "3"),
        "thermostat"   => ("Polynomial",     "1D Thermostat",        "3"),
        "oscillator"   => ("Polynomial",     "2D Oscillator",        "3"),
        "pendulum"     => ("PWA Inclusion",  "2D Pendulum",          "4"),
        "unicycle"     => ("PWA Inclusion",  "4D Unicycle",          "4"),
    )
    class_name, system_name, table = get(system_map, system_dir, (class_dir, system_dir, ""))

    barrier_type = barrier_dir == "PWC" ? "PWC" : "SOS"
    pwc_method = pwc_dir == "DUAL" ? "DUAL" :
                 pwc_dir == "CEGIS" ? "CEGIS" :
                 pwc_dir == "GD" ? "GD" : ""

    degree = ""
    if barrier_type == "SOS" && haskey(config, "barrier_settings") &&
       haskey(config["barrier_settings"], "barrier_degree")
        degree = string(config["barrier_settings"]["barrier_degree"])
    end

    # num_partitions: for linear PWC we can compute from ϵ × state_space, for
    # nonlinear (PWA inclusion) we read it from the yaml filename (SOS + PWC).
    num_partitions = ""
    if system_dir in ("pendulum", "unicycle")
        fname = splitext(basename(yaml_file))[1]
        m = match(r"(\d+)$", fname)
        num_partitions = m === nothing ? "" : m.captures[1]
    elseif barrier_type == "PWC" && system_dir == "contraction2"
        # partition count from epsilon × state_space
        try
            ϵ = config["transition_probalities"]["ϵ"]
            low  = config["state_space"]["low"]
            high = config["state_space"]["high"]
            n = prod(ceil.(Int, (high .- low) ./ ϵ))
            num_partitions = string(n)
        catch
            num_partitions = ""
        end
    elseif barrier_type == "SOS"
        num_partitions = "NA"
    end

    return (table, class_name, system_name, barrier_type, pwc_method, degree, num_partitions)
end

_fmt(::Missing) = ""
_fmt(x::AbstractString) = x
_fmt(x) = string(x)

function write_result_row(yaml_file::AbstractString, config::AbstractDict;
                          status::AbstractString,
                          time_s = "",
                          eta_val = "",
                          beta_val = "",
                          Ps_val = "")
    table, class_name, system_name, barrier_type, pwc_method, degree, num_partitions =
        parse_yaml_metadata(yaml_file, config)

    path = results_path()
    _ensure_header(path)
    open(path, "a") do io
        fields = ["julia", table, class_name, system_name, barrier_type, pwc_method,
                  degree, num_partitions,
                  _fmt(time_s), _fmt(eta_val), _fmt(beta_val), _fmt(Ps_val), status]
        println(io, join(fields, ","))
    end
    return nothing
end

function write_result_row(yaml_file::AbstractString, config::AbstractDict, res;
                          time_horizon::Integer)
    write_result_row(yaml_file, config;
        status = "OK",
        time_s = total_time(res),
        eta_val = eta(res),
        beta_val = beta(res),
        Ps_val = psafe(res, time_horizon))
end
