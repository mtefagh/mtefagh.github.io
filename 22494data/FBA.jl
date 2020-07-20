using JSON, SparseArrays
open("Recon3D.json") do f
    dict = JSON.parse(f);
    data = dict["reactions"];
    temp_name = [(data[i])["name"] for i in 1:length(data)];
    b = findall(temp_name .== "Generic Human Biomass Reaction")[1];
    order = [1:(b-1); (b+1):(length(temp_name)); b];
    global name = [(data[i])["name"] for i in order];
    global lower_bound = [(data[i])["lower_bound"] for i in order];
    global upper_bound = [(data[i])["upper_bound"] for i in order];
    global subsystem = [(data[i])["subsystem"] for i in order];
    metabolites = [(data[i])["metabolites"] for i in order];
    id = [((dict["metabolites"])[i])["id"] for i in 1:length(dict["metabolites"])];
    S = Matrix{Float64}(undef, length(id), length(metabolites));
    for i in 1:length(metabolites)
        for j in 1:length(id)
            S[j, i] = get(metabolites[i], id[j], 0)
        end
    end
    global S = sparse(S); 
end
