using Printf
include("./900_common_headers.jl")

struct TCM #TeraChem molden parsing tags
    outtag_m = "./"
    # Define start and end blocks for molden file
    coord_st = "Atoms"
    coord_ed = "GTO"
    gto_end  = "MO"
    mo_end   = ",\$p"
    root_str = "Ene"
    # Define output file names
    fcoord   = string(outtag_m,"coords.inp") 
    fgto     = string(outtag_m,"gto.inp") 
    fmovec   = string(outtag_m,"movec.inp") 
    froots   = string(outtag_m,"roots.inp") 
    febasis  = string(outtag_m,"ebasis.inp") 
    fmat     = string(outtag_m,"mo_to_ao.inp") 
    froot    = string(outtag_m,"moeng.inp") 
    focc     = string(outtag_m,"occupancy.inp") 
    fspin    = string(outtag_m,"spin.inp")
end # struct

function dp_format(args...)
    # Format args in double precision
    s = ""
    for a in args
        s *= @sprintf("%.15e ",a)
    end
    return s
end

function lnx_cmd(comm::Cmd)
    # Run linux shell command in script
    # The comm variable must be a comm object in julia, denoted by backticks '
    run(bash -c comm)
    return Nothing
end

function line_cnt_file(infile)
    # Get line number from a file
    lc = countlines(infile)
    return lc
end

function grep_str_to_file(string,fin,fout)
    string = "$string"
    open(fout,"w") do io
        run(pipeline(`grep $string $fin`,stdout=io))
    end
    close(fout)
    return Nothing
end

function extract_lines_from_file(pos_start, pos_end, infile, outfile, EOF_ext=true)
    ftmp  = "./tmp.sed"
    ftmp2 = "./tmp.cat"
    comm1 = `sed -n $("/" * pos_start * "/,/" * pos_end * "/p") $infile`
    comm2 = `sed "\$d" $ftmp`
    comm3 = `tail -n+2 $ftmp2`
    if EOF_ext == true
        run(pipeline(comm1,stdout=ftmp))
        run(pipeline(comm2,stdout=ftmp2))
        run(pipeline(comm3,stdout=fout))
        rm(ftmp2)
    else
        comm1 = `sed -n $("/" * pos_start * "/" * pos_end) $infile`
        comm3 = `tail -n+2 $ftmp`
        run(pipeline(comm1,stdout=ftmp))
        run(pipeline(comm3,stdout=fout))
        rm(ftmp)
    end
    return Nothing
end

function extract_momat_molden(nb,nao,infile)
    mat   = zeros((nb,nb))
    roots = zeros(nb)
    fin   = open(infile,"r")
    for imo=1:nb
        ene_line   = readline(fin)
        ene_line   = strip(split(ene_line))
        roots[imo] = parse(Float64,ene_line[1])
        spin_line  = readline(fin)
        occup_line = readline(fin)
        for iao=1:nao
            line = readline(fin)
            line = strip(split(line))
            mat[iao,imo] = parse(Float64,line[1])
        end
    end
    close(fin)
    return roots,mat
end

function read_last_mo(imo,nb,infile)
    vec = zeros(nb)
    open(infile,"r")
    for i=1:nb
        line = readline(infile)
        line = strip(split(line))
        if length(line) == 2
            iao = parse(Int,line[1]) - 1
            vec[iao] = parse(Float64,line[1])
        else
            break
        end
    end
    close(infile)
    return vec
end













