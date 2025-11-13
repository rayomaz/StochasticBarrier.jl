using NCDatasets

# === Specify the filename ===
# filename = "pendulum/probability_data_120_sigma_[0.1, 0.1].nc"
filename = "unicycle/probability_data_1250_sigma_[0.01, 0.01, 0.01, 0.01].nc"

# === Open dataset explicitly from NCDatasets ===
ds = NCDatasets.Dataset(filename, "r")

# === List all variable names ===
println("Variables in file:")
varnames = collect(keys(ds))
println(varnames)

# === Read all variables into a dictionary ===
data = Dict{String, Any}()

for varname in varnames
    println("Reading variable: ", varname)
    data[string(varname)] = ds[varname][:]  # load entire variable
end

# === Close file ===
close(ds)

# === Display summary of contents ===
println("\nSummary of loaded variables:")
for (name, val) in data
    println(rpad(name, 35), ": size = ", size(val), ", eltype = ", eltype(val))
end
