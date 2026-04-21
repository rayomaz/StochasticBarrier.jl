#!/usr/bin/env julia
# Parse per-tool CSV outputs into the two summary tables from the paper.
#
#   $SB_RESULTS_DIR/julia_results.csv     (Julia)
#   $SB_RESULTS_DIR/protect_results.csv   (PRoTECT)
#   $SB_RESULTS_DIR/sostools_results.csv  (SOSTOOLS)
#
# writes
#   $SB_RESULTS_DIR/table3.csv
#   $SB_RESULTS_DIR/table4.csv
#
# Usage: julia --project=. parse_results.jl [results_dir]

module ParseResults

using Printf

const COLS = (:tool, :table, :class, :system, :barrier_type, :pwc_method,
              :degree, :num_partitions, :time_s, :eta, :beta, :Ps, :status)

struct Row
    tool::String
    table::String
    class::String
    system::String
    barrier_type::String
    pwc_method::String
    degree::String
    num_partitions::String
    time_s::String
    eta::String
    beta::String
    Ps::String
    status::String
end

function _split_csv_line(line::AbstractString)
    # Our writers never quote fields (numbers + simple words), so a split is fine.
    return String.(split(line, ','))
end

function read_csv(path::AbstractString)::Vector{Row}
    rows = Row[]
    isfile(path) || return rows
    open(path, "r") do io
        header_line = readline(io)
        header = _split_csv_line(header_line)
        idx = Dict(Symbol(strip(h)) => i for (i, h) in enumerate(header))
        for line in eachline(io)
            isempty(strip(line)) && continue
            parts = _split_csv_line(line)
            length(parts) == length(header) || continue
            get_col(k) = get(idx, k, 0) == 0 ? "" : strip(parts[idx[k]])
            push!(rows, Row(
                get_col(:tool), get_col(:table), get_col(:class), get_col(:system),
                get_col(:barrier_type), get_col(:pwc_method),
                get_col(:degree), get_col(:num_partitions),
                get_col(:time_s), get_col(:eta), get_col(:beta), get_col(:Ps),
                get_col(:status)))
        end
    end
    return rows
end

# -------------- formatting ----------------

function fmt_time(s::AbstractString, status::AbstractString)
    status == "OM" && return "OM"
    status == "FAILED" && return tryparse(Float64, s) === nothing ? "FAILED" :
                                 @sprintf("%.2f", parse(Float64, s))
    v = tryparse(Float64, s)
    v === nothing && return s
    return @sprintf("%.2f", v)
end

function fmt_sci(s::AbstractString, status::AbstractString; default="×")
    status == "OM" && return "×"
    status == "FAILED" && return "FAILED"
    v = tryparse(Float64, s)
    v === nothing && return default
    isnan(v) && return default
    if v == 0
        return "0.00"
    end
    if abs(v) >= 0.1
        return @sprintf("%.2f", v)
    end
    # 1.7e-2 style
    exp10v = floor(Int, log10(abs(v)))
    mant = v / 10.0^exp10v
    return @sprintf("%.1fe%d", mant, exp10v)
end

function fmt_Ps(s::AbstractString, status::AbstractString)
    status == "OM" && return "×"
    status == "FAILED" && return "0.00"
    v = tryparse(Float64, s)
    v === nothing && return "×"
    return @sprintf("%.2f", v)
end

function _find(rows::Vector{Row}; tool, system, degree=nothing,
               barrier_type=nothing, pwc_method=nothing, num_partitions=nothing)
    for r in rows
        r.tool == tool || continue
        r.system == system || continue
        degree === nothing || r.degree == string(degree) || continue
        barrier_type === nothing || r.barrier_type == barrier_type || continue
        pwc_method === nothing || r.pwc_method == pwc_method || continue
        num_partitions === nothing || r.num_partitions == string(num_partitions) || continue
        return r
    end
    return nothing
end

# ---- Table 3 ----

const TABLE3_ROWS = [
    # (class, system, [degrees])
    ("Linear",     "2D Contraction Map 1", [2, 4, 8, 12, 24, 30]),
    ("Linear",     "2D Two Tank",          [10, 12, 14]),
    ("Linear",     "3D Room Temp.",        [6, 8, 14]),
    ("Linear",     "6D Quadrotor",         [2, 4, 6]),
    ("Polynomial", "1D Thermostat",        [6, 8, 10]),
    ("Polynomial", "2D Oscillator",        [8, 12, 14]),
]

function build_table3(rows::Vector{Row})
    header = ["Class", "System", "Deg.",
              "sostools_tau", "sostools_eta", "sostools_beta", "sostools_Ps",
              "protect_tau",  "protect_eta",  "protect_beta",  "protect_Ps",
              "julia_tau",    "julia_eta",    "julia_beta",    "julia_Ps"]
    out = String[join(header, ",")]

    for (cls, sys, degs) in TABLE3_ROWS
        for (i, d) in enumerate(degs)
            cls_col = i == 1 ? cls : ""
            sys_col = i == 1 ? sys : ""

            cells = [cls_col, sys_col, string(d)]
            for tool in ("sostools", "protect", "julia")
                r = _find(rows; tool=tool, system=sys, degree=d, barrier_type="SOS")
                status = r === nothing ? "OM" : r.status
                push!(cells, fmt_time(r === nothing ? "" : r.time_s, status))
                push!(cells, fmt_sci(r === nothing ? "" : r.eta, status))
                push!(cells, fmt_sci(r === nothing ? "" : r.beta, status))
                push!(cells, fmt_Ps(r === nothing ? "" : r.Ps, status))
            end
            push!(out, join(cells, ","))
        end
    end
    return join(out, "\n") * "\n"
end

# ---- Table 4 ----
# All-Julia: SOS vs PWC (DUAL, CEGIS, GD).

const TABLE4_ROWS = [
    # (class, system, [(sos_Q, sos_deg, pwc_Q, pwc_method_Q_override)])
    # For Contraction2 the SOS |Q| is "NA" and PWC |Q| differs per row.
    ("Linear",        "2D Contraction Map 2", [
        (nothing, 4,  nothing),   # PWC partition count filled below from rows
        (nothing, 8,  nothing),
        (nothing, 20, nothing),
        (nothing, 30, nothing),
    ]),
    ("PWA Inclusion", "2D Pendulum", [
        (120, 2, 120),
        (240, 2, 240),
        (480, 2, 480),
    ]),
    ("PWA Inclusion", "4D Unicycle", [
        (1250, 4, 1250),
        (1800, 4, 1800),
    ]),
]

function _pwc_contraction2_partitions(rows::Vector{Row}, method::String)
    # Contraction2 PWC rows have varying |Q|. We pair degree-sorted with partition-sorted.
    # Return vector of (num_partitions::String, Row) sorted by num_partitions ascending.
    pairs = Tuple{Int,Row}[]
    for r in rows
        r.tool == "julia" || continue
        r.system == "2D Contraction Map 2" || continue
        r.barrier_type == "PWC" || continue
        r.pwc_method == method || continue
        n = tryparse(Int, r.num_partitions)
        n === nothing && continue
        push!(pairs, (n, r))
    end
    sort!(pairs; by = first)
    return pairs
end

function build_table4(rows::Vector{Row})
    header = ["Class", "System",
              "SOS_Q", "SOS_Deg",  "SOS_tau", "SOS_Ps",
              "DUAL_Q", "DUAL_tau", "DUAL_Ps",
              "CEGIS_tau", "CEGIS_Ps",
              "GD_tau",    "GD_Ps"]
    out = String[join(header, ",")]

    # Pre-compute contraction2 PWC partition ordering per method
    c2_dual  = _pwc_contraction2_partitions(rows, "DUAL")
    c2_cegis = _pwc_contraction2_partitions(rows, "CEGIS")
    c2_gd    = _pwc_contraction2_partitions(rows, "GD")

    for (cls, sys, entries) in TABLE4_ROWS
        for (i, (sos_Q, sos_deg, pwc_Q)) in enumerate(entries)
            cls_col = i == 1 ? cls : ""
            sys_col = i == 1 ? sys : ""

            # SOS
            rS = sos_Q === nothing ?
                _find(rows; tool="julia", system=sys, barrier_type="SOS", degree=sos_deg) :
                _find(rows; tool="julia", system=sys, barrier_type="SOS",
                      degree=sos_deg, num_partitions=sos_Q)
            sos_status = rS === nothing ? "OM" : rS.status
            sos_Q_str = sos_Q === nothing ? "NA" : string(sos_Q)

            # PWC — for contraction2, pick the i-th entry per method.
            function pwc_cell(method_rows, explicit_Q)
                if sys == "2D Contraction Map 2"
                    idx = i
                    if idx ≤ length(method_rows)
                        n, r = method_rows[idx]
                        return (string(n), r)
                    end
                    return ("", nothing)
                else
                    r = _find(rows; tool="julia", system=sys, barrier_type="PWC",
                              num_partitions=explicit_Q)
                    Qstr = explicit_Q === nothing ? "" : string(explicit_Q)
                    return (Qstr, r)
                end
            end

            dual_Q, rD = sys == "2D Contraction Map 2" ? pwc_cell(c2_dual, nothing) :
                                                         ("", _find(rows; tool="julia", system=sys,
                                                                     barrier_type="PWC",
                                                                     pwc_method="DUAL",
                                                                     num_partitions=pwc_Q))
            _, rC = sys == "2D Contraction Map 2" ? pwc_cell(c2_cegis, nothing) :
                                                    ("", _find(rows; tool="julia", system=sys,
                                                                barrier_type="PWC",
                                                                pwc_method="CEGIS",
                                                                num_partitions=pwc_Q))
            _, rG = sys == "2D Contraction Map 2" ? pwc_cell(c2_gd, nothing) :
                                                    ("", _find(rows; tool="julia", system=sys,
                                                                barrier_type="PWC",
                                                                pwc_method="GD",
                                                                num_partitions=pwc_Q))

            # For non-contraction2 the |Q| cell matches the pwc_Q hint.
            if sys != "2D Contraction Map 2"
                dual_Q = pwc_Q === nothing ? "" : string(pwc_Q)
            end

            cells = [cls_col, sys_col,
                     sos_Q_str, string(sos_deg),
                     fmt_time(rS === nothing ? "" : rS.time_s, sos_status),
                     fmt_Ps(rS === nothing ? "" : rS.Ps, sos_status)]

            # DUAL (with leading |Q|)
            status = rD === nothing ? "OM" : rD.status
            push!(cells, dual_Q)
            push!(cells, fmt_time(rD === nothing ? "" : rD.time_s, status))
            push!(cells, fmt_Ps(rD === nothing ? "" : rD.Ps, status))
            # CEGIS
            status = rC === nothing ? "OM" : rC.status
            push!(cells, fmt_time(rC === nothing ? "" : rC.time_s, status))
            push!(cells, fmt_Ps(rC === nothing ? "" : rC.Ps, status))
            # GD
            status = rG === nothing ? "OM" : rG.status
            push!(cells, fmt_time(rG === nothing ? "" : rG.time_s, status))
            push!(cells, fmt_Ps(rG === nothing ? "" : rG.Ps, status))

            push!(out, join(cells, ","))
        end
    end
    return join(out, "\n") * "\n"
end

function main(args::Vector{String} = ARGS)
    dir = isempty(args) ? get(ENV, "SB_RESULTS_DIR",
                              joinpath(@__DIR__, "results")) : args[1]
    isdir(dir) || error("results directory not found: $dir")

    rows = vcat(
        read_csv(joinpath(dir, "julia_results.csv")),
        read_csv(joinpath(dir, "protect_results.csv")),
        read_csv(joinpath(dir, "sostools_results.csv")),
    )

    t3 = joinpath(dir, "table3.csv")
    t4 = joinpath(dir, "table4.csv")
    write(t3, build_table3(rows))
    write(t4, build_table4(rows))
    @info "Wrote tables" table3=t3 table4=t4 rows=length(rows)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    ParseResults.main(ARGS)
end
