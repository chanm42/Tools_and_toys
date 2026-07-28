using Printf
include("./900_common_header.jl")

struct TCM
    outtag_m::String
    coord_st::String
    coord_ed::String
    gto_end ::String
    mo_end  ::String
    root_str::String
    fcoord  ::String
    fgto    ::String
    fmovec  ::String
    froots  ::String
    febasis ::String
    fmat    ::String
    froot   ::String
    focc    ::String
    fspin   ::String
end # struct

function TCM_cons(outtag_m) #TeraChem molden parsing tag
    # Define start and end blocks for molden file
    coord_st = "Atoms"
    coord_ed = "GTO"
    gto_end  = "MO"
    mo_end   = ",\$p"
    root_str = "Ene"
    # Define output file names
    fcoord   = joinpath(outtag_m,"coords.inp") 
    fgto     = joinpath(outtag_m,"gto.inp") 
    fmovec   = joinpath(outtag_m,"movec.inp") 
    froots   = joinpath(outtag_m,"roots.inp") 
    febasis  = joinpath(outtag_m,"ebasis.inp") 
    fmat     = joinpath(outtag_m,"mo_to_ao.inp") 
    froot    = joinpath(outtag_m,"moeng.inp") 
    focc     = joinpath(outtag_m,"occupancy.inp") 
    fspin    = joinpath(outtag_m,"spin.inp")
    obj      = TCM(outtag_m,
                   coord_st,
                   coord_ed,
                   gto_end,
                   mo_end,
                   root_str,
                   fcoord,
                   fgto,
                   fmovec,
                   froots,
                   febasis,
                   fmat,
                   froot,
                   focc,
                   fspin)
    return obj
end # constructor

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
    # Typically easier to use the run(pipeline(comm,stdout=XXX)) method if piping
    run(`bash -c $comm`)
    return nothing
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
    return nothing
end

function extract_lines_from_file(pos_start::String, pos_end::String, infile::String, outfile::String, EOF_ext=true)
    ftmp  = "./tmp.sed"
    ftmp2 = "./tmp.cat"
    comm1 = `sed -n $("/" * pos_start * "/,/" * pos_end * "/p") $infile`
    comm2 = `sed "\$d" $ftmp`
    comm3 = `tail -n+2 $ftmp2`
    if EOF_ext == true
        run(pipeline(comm1,stdout=ftmp))
        run(pipeline(comm2,stdout=ftmp2))
        run(pipeline(comm3,stdout=outfile))
        rm(ftmp2)
    else
        comm1 = `sed -n $("/" * pos_start * "/" * pos_end) $infile`
        comm3 = `tail -n+2 $ftmp`
        run(pipeline(comm1,stdout=ftmp))
        run(pipeline(comm3,stdout=outfile))
        rm(ftmp)
    end
    return nothing
end

function extract_momat_v2(nb,nao,infile)
    lc    = countlines(infile)
    mat   = zeros((nb,nb))
    roots = zeros(nb)
    spin  = zeros(nb)
    occ   = zeros(nb)
    fin   = open(infile,"r")
    imo::Int = 0 
    for l=1:lc
        line = readline(fin)
        if occursin("ENE",uppercase(line))
            imo += 1
            if imo != nb
                e = split(strip(line))
                roots[imo] = parse(Float64,e[2])
                continue
            end
        elseif occursin("SPIN",uppercase(line))
            if occursin("ALPHA",uppercase(line))
                spin[imo] = 0.50e0
            elseif occursin("BETA",uppercase(line))
                spin[imo] = -0.50e0
            end
            continue
        elseif occursin("OCCUP",uppercase(line))
            oc = split(strip(line))
            occ[imo] = parse(Float64,oc[2])
            continue
        elseif imo == nb
            break
        elseif occursin("MOLDEN",uppercase(line))
            break
        else
            lin  = split(strip(line))
            iao  = parse(Int,lin[1])
            mval = parse(Float64,lin[2])
            mat[iao,imo] = mval
        end
    end
    close(fin)
    return spin,occ,roots,mat
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

function write_mat(nb,mat,filename)
    fout  = open(filename,"w")
    buff1 = "NAO \t NMO \n"
    buff2 = "$nb \t $nb \n"
    buff3 = "IAO \t IMO \t COEFF[IAO,IMO]\n"
    write(fout,buff1)
    write(fout,buff2)
    write(fout,buff3)

    for imo=1:nb
        for iao=1:nb
            line = "$iao \t $imo \t $(dp_format(mat[iao,imo])) \n"
            write(fout,line)
        end
    end
    close(fout)
    return nothing
end

function write_vec(nb,vec,filename)
    fout = open(filename,"w")
    line = "$nb \n"
    write(fout,line)
    for imo=1:nb
        line = "$imo \t $(dp_format(vec[imo]))\n"
        write(fout,line)
    end
    close(fout)
    return nothing
end

function gto_get_atom_breaks(fgto)
    comm1 = `sed -n '/^$/=' $fgto`
    run(pipeline(comm1,stdout="tmp.break"))
    lc    = countline(tmp.break)
    atm_start  = 0
    atm_breaks = []
    fin   = open(tmp.break,"r")
    for idx=1:lc
        line = readline(tmp.break)
        line = split(strip(line))
        atm_break = [atm_start,parse(Int,line[1])]
        push!(atm_breaks,atm_break)
        atm_start = parse(Int,line[1])
    end
    rm(tmp.break)
    
    return atm_breaks
end

function read_one_bas!(lines_iter,lsym,ncont,fac=1.0)
    fac = parse(Float64,fac)
    bas = Vector{Any}[]
    nc  = parse(Int,ncont)
    for i=1:nc
        line = popfirst!(lines_iter)
        dat  = split(line)
        push!(bas,
              [
               uppercase(String(lsym)),
               nc,
               i-1,
               parse(Float64,dat[1]),
               parse(Float64,dat[2])*fac
              ])
    end
    return bas
end

function molden_gto_reader(natm,coords,fgto)
    fin        = open(fgto,"r")
    lines      = readlines(fin)
    close(fin)
    lines_iter = Iterators.Stateful(lines)
    cent       = Vector{Any}[]
    basis      = Dict{Int,Any}()
    iatm::Int  = 0
    
    for line in lines_iter
        dat = split(line)
        if isempty(dat)
            continue
        elseif tryparse(Int,dat[1]) !== nothing
            #println(">>Tryparse sucsessfull= $(dat[1])")
            iatm  = parse(Int,dat[1])
            coord = coords[iatm]
            push!(cent,[iatm,coord])
            basis[iatm] = []
            #println("cent= $cent")
        elseif uppercase(dat[1]) in ["S","P","D","F","G","H","I","J"]
            #println("dat= $dat")
            push!(basis[iatm],read_one_bas!(
                                           lines_iter,
                                           dat[1],
                                           dat[2],
                                           dat[3]))
        end # if
    end # lines
    return cent,basis
end

function ebasis_writer(natm,centl,basis,fileout)
    fout = open(fileout,"w")
    line = "NO. OF UNIQUE CENTERS \n"
    write(fout,line)
    line = "$natm \n"
    write(fout,line)
    for iatm=1:natm
        icent  = centl[iatm][2]
        buff   = "---------------$(icent[2])--------------- \n"
        write(fout,buff)
        icoord = icent[3]
        if false
            println("icent= $icent")
            println("length= $(length(icent))")
            println("icoord= $icoord")
            println("length= $(length(icoord))")
            println("centl= $centl")
        end
        s      = dp_format(icoord[1],icoord[2],icoord[3])
        buff1  = "CENTER \n"
        buff2  = "$s \n"
        buff3  = "SYMBOLS \n"
        buff5  = "  $(string(length(basis[iatm]))) \n"
        write(fout,buff1)
        write(fout,buff2)
        write(fout,buff3)
        write(fout,buff5)
        for iprim=1:length(basis[iatm])
            sym   = basis[iatm][iprim][1][1]
            ncont = basis[iatm][iprim][1][2]
            buff  = "$sym \t $(string(ncont)) \n"
            write(fout,buff)
            for icont=1:ncont
                alpha = basis[iatm][iprim][icont][4]
                coeff = basis[iatm][iprim][icont][5]
                #println("raw_alpha= $alpha")
                #println("raw_coeff= $coeff")
                #println("dp_alpha= $(dp_format(alpha))")
                #println("dp_coeff= $(dp_format(coeff))")
                buff  = "\t $(string(icont)) \t $(dp_format(alpha)) \t $(dp_format(coeff)) \n"
                #println(buff)
                write(fout,buff)
            end # for icont
        end # for iprim
    end # for iatm
    close(fout)
    return nothing
end

function read_coords(natm,fcoord)
    ang_to_bohr = 1.8897259886e0
    fin    = open(fcoord,"r")
    coords = Vector{Any}[]
    for i=1:natm
        line = readline(fin)
        line = split(strip(line))
        atm  = string(line[1])
        iatm = parse(Int,line[2])
        inuc = parse(Int,line[3])
        x    = parse(Float64,line[4])*ang_to_bohr
        y    = parse(Float64,line[5])*ang_to_bohr
        z    = parse(Float64,line[6])*ang_to_bohr
        xyz  = [x,y,z]
        push!(coords,xyz)
    end
    close(fin)
    return coords
end

function read_coords_v2(natm,fcoord)
    ang_to_bohr = 1.8897259886e0
    fin    = open(fcoord,"r")
    coords = Vector{Any}[]
    xmin   = zeros(3)
    xmax   = zeros(3)
    for i=1:3
        xmin[i] = 1.0e6
        xmax[i] = -1.0e6
    end
    for i=1:natm
        line = readline(fin)
        line = split(strip(line))
        atm  = string(line[1])
        iatm = parse(Int,line[2])
        inuc = parse(Int,line[3])
        x    = parse(Float64,line[4])*ang_to_bohr
        y    = parse(Float64,line[5])*ang_to_bohr
        z    = parse(Float64,line[6])*ang_to_bohr
        xyz  = [x,y,z]
        for idim=1:3
            if xmin[idim] > xyz[idim]
                xmin[idim] = xyz[idim]
            elseif xmax[idim] < xyz[idim]
                xmax[idim] = xyz[idim]
            end
        end # for idim
        push!(coords,[iatm,atm,xyz])
    end # iatm
    close(fin)
    return xmin,xmax,coords
end

function test_momat(nb,roots,momat)
    fout = open("test_momat.inp","w")
    for i=1:nb
        str1 = "Ene= $(dp_format(roots[i])) \n"
        write(fout,buff)
        str2 = "SPIN \n"
        str3 = "OCCUP \n"
        write(fout,str2)
        write(fout,str3)
        for j=1:nb
            str = "\t $(string(j)) $(dp_format(momat[i,j])) \n"
            write(fout,str)
        end
    end
    close(fout)
    return nothing
end

function parse_mo(filename,outtag)
    println("---- Parsing $filename ----")
    obj = TCM_cons(outtag)

    extract_lines_from_file(obj.coord_st,obj.coord_ed,filename,obj.fcoord)
    extract_lines_from_file(obj.coord_ed,obj.gto_end,filename,obj.fgto)
    extract_lines_from_file(obj.gto_end,obj.mo_end,filename,obj.fmovec,false)
    grep_str_to_file(obj.root_str,obj.fmovec,obj.froots)
    nb = countlines(obj.froots)
    spin,occ,roots,momat = extract_momat_v2(nb,nb,obj.fmovec)
    write_mat(nb,momat,obj.fmat)
    write_vec(nb,roots,obj.froot)
    write_vec(nb,spin,obj.fspin)
    write_vec(nb,occ,obj.focc)

    natm = countlines(obj.fcoord)
    xmin,xmax,coords = read_coords_v2(natm,obj.fcoord)
    centl, basis     = molden_gto_reader(natm,coords,obj.fgto)
    ebasis_writer(natm,centl,basis,obj.febasis)

    return xmin,xmax,nb
end

function write_grid_maker(xmin,xmax,nb,filename)
    fout = open(filename,"w")
    str1 = "NDIM \n"
    str2 = "3 \n"
    write(fout,str1)
    write(fout,str2)
    str3 = "IDIM NPT_IDIM IMIN \t IMAX \n"
    write(fout,str3)
    for i=1:3
        str = "$i \t 100 \t $(dp_format(xmin[i])) \t $(dp_format(xmax[i])) \n"
        write(fout,str)
    end
    str4 = "NMO_for_grid NMO_MAX \n"
    str5 = "XXX \t $nb \n"
    str6 = "IDX IMO \n"
    write(fout,str4)
    write(fout,str5)
    write(fout,str6)
    close(fout)
    return nothing
end

function test1A()
    # filename = path to molden file
    # outtag   = directory you want files published to
    filename = "./b10n10s.molden"
    outtag   = "./"
    xmin,xmax,nb = parse_mo(filename,outtag)
    println("Finished parsing $filename")
    println("> xmax= $xmax")
    println("> xmin= $xmin")
    delta = 12.0e0
    xmin  = xmin .- delta
    xmax  = xmax .+ delta
    fgrid = outtag*"grid_maker_info.inp"
    write_grid_maker(xmin,xmax,nb,fgrid)
    
    return nothing
end


if abspath(PROGRAM_FILE) == @__FILE__
  #---
    println()
    println("Starting main from file: ", PROGRAM_FILE)
    println("--------------------------------------------")
    println("start_time= ",now())
    println("--------------------------------------------")
    println()
  #---
        test1A()
  #---
    println()
    println("--------------------------------------------")
    println("end_time= ",now())
    println("--------------------------------------------")
  #---
end #abspath

