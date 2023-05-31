using ArgParse

function parse_commandline(args_array::Array{String,1}, appfolder::String)
    s = ArgParseSettings(usage="##### Convert arguments #####",exit_after_help=false)
    @add_arg_table s begin
       "instance"
        help = "Instance file path"
        "--out","-o"
         help = "Path to write the instance"
        "--lvalue","-l" #time limit for each trip
         help = "l value - new value for the time limit for each trip"
         arg_type = Int
    end
    return parse_args(args_array, s)
 end


 function writeinstance(solpath, aux)
    open(solpath, "w") do f
        write(f, aux[1], " ", aux[2], " ", aux[3], "\n") 
        pos = 4
        for i=1:parse(Int, aux[1])
            write(f, aux[pos], " ", aux[pos+1], " ", aux[pos+2],"\n") 
            pos +=3
        end
        for i=1:parse(Int, aux[2])
            write(f, aux[pos], " ", aux[pos+1], " ", aux[pos+2]," ", aux[pos+3],"\n") 
            pos +=4
        end

    end
 end

function main()
    appfolder = dirname(@__FILE__)
    app = parse_commandline(ARGS, appfolder)
    str = read(app["instance"], String)
    breaks_in = [' '; ':'; '\n';'\t';'\r']
    aux = split(str, breaks_in; limit=0, keepempty=false)
    aux[3] = string(app["lvalue"])
    writeinstance(app["out"],aux)
end

main()