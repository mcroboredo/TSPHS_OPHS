
# Leitura dos valores ótimos e criação de um dicionário chamado otimos_dict
file = open("otimos.txt", "r")
dados = read(file, String)
breaks_in = [' '; ':'; '\n';'\t';'\r']
aux = split(dados, breaks_in; limit=0, keepempty=false)    
otimos_dict = Dict()
nrows = Int(length(aux)/2)
for i = 1:2:nrows
    push!(otimos_dict, aux[i] => aux[i+1]) 
end 
close(file)

# Leitura do tour ótimo
file = open("Data_tours/eil76.opt.tour", "r")
dados = read(file, String)
breaks_in = [' '; ':'; '\n';'\t';'\r']
aux = split(dados, breaks_in; limit=0, keepempty=false)
nrows = length(aux)
sol = []
for i = 1:nrows
    if aux[i] == "TOUR_SECTION"
        for j = i+1:length(aux)-1
            push!(sol, parse(Int,aux[j]))
        end
    end
end
close(file)

# leitura do arquivo da instância

