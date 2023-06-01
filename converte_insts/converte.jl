lista_pastas = readdir("Dados_para_converter")
for pasta in lista_pastas
    lista_arquivos = readdir("Dados_para_converter/" * pasta)
    for arquivo in lista_arquivos
        caminho_arquivo = "Dados_para_converter/" * pasta  * "/" * arquivo 
        # Read the file and obtain aux
        file = open(caminho_arquivo, "r")
        dados = read(file, String)
        close(file)
        breaks_in = [' '; ':'; '\n';'\t';'\r']
        aux = split(dados, breaks_in; limit=0, keepempty=false)
        # Write other file
        caminho_novo = "Dados_convertidos/" * pasta  * "/" * arquivo
        open(caminho_novo, "w") do f
            n_hotels = parse(Int, aux[2]) + 2 
            n_customers = parse(Int, aux[1]) - 2
            write(f, string(n_hotels), " ", string(n_customers), " ", aux[5], "\n") # Write first row (|H| |C| |T_d|)
            # Loop to write the hotels
            pos = 4 + parse(Int,aux[3]) + 1 # Get the position of the first hotel
            count = 0
            for i=1:n_hotels
                write(f, string(count), " ", aux[pos], " ", aux[pos+1], "\n")
                pos += 3
                count += 1
            end
            # Loop to write clients info
            for i=1:n_customers
                write(f, string(count), " ", aux[pos], " ", aux[pos+1], " ", "0", "\n")
                pos += 3
                count += 1
            end
        end
    end
end

