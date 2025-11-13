using NCDatasets

# filename = "pendulum/dynamics_120.nc"
filename = "unicycle/dynamics_1250.nc"

# Open the dataset (read-only)
ds = NCDatasets.Dataset(filename, "r")

# === Display all variable names ===
println("Variables in file:")
varnames = keys(ds)
println(collect(varnames))

# === Read each variable into a dictionary ===
data = Dict{String, Any}()

for varname in varnames
    println("Reading variable: ", varname)
    data[string(varname)] = ds[varname][:]   # read all data
end

# Close the dataset
close(ds)

# === Inspect what's loaded ===
for (name, val) in data
    println(rpad(name, 25), ": size = ", size(val), ", eltype = ", eltype(val))
end
