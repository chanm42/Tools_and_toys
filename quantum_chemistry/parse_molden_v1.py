import numpy as np
import scipy
import datetime
import sys
import subprocess
import scipy.linalg as la

class TCM(object): #TeraChem Molden Parsing tags
    def __init__(self,outdir='./'):
        self.coord_st = 'Atoms'
        self.coord_ed = 'GTO'
        self.gto_end  = 'MO'
        self.mo_end   = ',$p'
        self.root_str = 'Ene'
        self.outdir   = outdir
        self.fcoord   = self.outdir+ 'coords.inp'
        self.fgto     = self.outdir+ 'gto.inp'
        self.fmovec   = self.outdir+ 'movec.inp'
        self.froots   = self.outdir+ 'roots.inp'
        self.febasis  = self.outdir+ 'ebasis.inp'
        self.fmat     = self.outdir+ 'mo_to_ao.inp'
        self.froot    = self.outdir+ 'moeng.inp'
        self.focc     = self.outdir+ 'occupancy.inp'
        self.fspin    = self.outdir+ 'spin.inp'

def dp_format(*args):
    s = ''
    for a in args:
        s = s + format(a,".15e") + " "
    return s

def lnx_cmd(comm):
    subprocess.call(comm, shell=True)
    return None

def line_cnt_file(filein):
    comm1 = 'wc -l '+filein + " > ./TMP.wc"
    lnx_cmd(comm1)
    fout = open('TMP.wc','r')
    line = fout.readline()
    line = line.strip().split()
    lc = int(line[0])
    fout.close()
    remove_tmp_file('TMP.wc')

    return lc

def grep_str_to_file(string,fin,fout):
    comm = 'grep "'+string+'" '+fin+' > ./'+fout
    lnx_cmd(comm)
    return None

def remove_tmp_file(ftmp):
    comm1 = "rm ./" + ftmp
    lnx_cmd(comm1)
    return None

def extract_lines_from_file(pos_start, pos_end, filename,fout,EOF_ext=True):
    ftmp = "TMP.sed"
    ftmp2 = "TMP.cat"
    comm1 = "sed -n " + '\'' + '/' + pos_start + '/,/' + pos_end + '/p' + '\'' + " " + filename + " > ./" + ftmp
    comm2 = "sed \\$d " + "./" + ftmp + " > ./" + ftmp2
    comm3 = "tail -n+2 ./" +ftmp2 + " > ./" + fout
    if EOF_ext == True:
        lnx_cmd(comm1)
        lnx_cmd(comm2)
        lnx_cmd(comm3)
        remove_tmp_file(ftmp2)
    else:
        comm1 = "sed -n " + '\'' + '/' + pos_start + '/' + pos_end +   '\'' + " " + filename + " > ./" + ftmp
        comm3 = "tail -n+2 ./" +ftmp + " > ./" + fout
        lnx_cmd(comm1)
        lnx_cmd(comm3)
        remove_tmp_file(ftmp)
    return None


def extract_momat_molden(nb,nao,filename):

    mat = np.zeros((nb,nb))
    roots = np.zeros(nb)
    fin = open(filename,'r')
    for imo in range(nb):
        ene_line   = fin.readline()
        ene_line = ene_line.strip().split() 
        roots[imo] = float(ene_line[1])
        spin_line  = fin.readline()
        occup_line = fin.readline()
        for iao in range(nao):
            line = fin.readline()
            line = line.strip().split()
            mat[iao,imo] = float(line[1])
    fin.close()

    return roots,mat

def read_last_mo(imo,nb,fin):
    vec = np.zeros(nb)

    for i in range(nb):
        line = fin.readline()
        line = line.strip().split()
        if len(line) == 2:
            iao  = int(line[0])-1
            vec[iao] = float(line[1])
        if len(line) != 2:
            break
    return vec

def extract_momat_v2(nb,nao,filename):

    fin = open(filename,'r')
    mat   = np.zeros((nb,nb))
    roots = np.zeros(nb)
    spin  = np.zeros(nb)
    occ   = np.zeros(nb)
    imo = -1
    for line in fin:
        if 'ENE' in line.upper():
            imo += 1
            if imo == nb:
                break
            else:
                roots[imo] = float(line.strip().split()[1])
                continue
        if 'SPIN' in line.upper():
            if 'ALPHA' in line.upper():
                spin[imo] = 0.50e0
            if 'BETA' in line.upper():
                spin[imo] = -0.50e0
            continue
        if 'OCCUP' in line.upper():
            occ[imo] = float(line.strip().split()[1])
            continue
        if 'MOLDEN' in line.upper():
            break
        else:
            iao = int(line.strip().split()[0]) - 1
            mat[iao,imo] = float(line.strip().split()[1])
            
    return spin,occ,roots,mat

def write_mat(nb,mat,filename):

    fout = open(filename,'w')
    buff = 'NAO    NMO' + '\n'
    buff += str(nb) + '\t'+str(nb)+ '\n'
    #buff = str(nb) + '\n'
    buff += "IAO   IMO   COEFF[IAO,IMO]" + '\n'
    fout.write(buff)

    for imo in range(nb):
        for iao in range(nb):
            buff = str(iao+1) +'\t'+ str(imo+1) + '\t' + dp_format(mat[iao,imo]) + '\n'
            fout.write(buff)
    fout.close()
    return None

def write_vec(nb,vec,filename):
    fout = open(filename,'w')
    #buff = 'NMO' + '\n'
    buff = str(nb) + '\n'
    #buff += "IMO   ENG[IMO]" + '\n'
    fout.write(buff)

    for imo in range(nb):
        buff = str(imo+1) + '\t' + dp_format(vec[imo]) + '\n'
        fout.write(buff)
    fout.close()

    return None

def gto_get_atom_breaks(fgto):
    
    comm = "sed -n '/^$/=' " + fgto + " > tmp.break"
    lnx_cmd(comm)
    fin = open('tmp.break','r')
    lines = fin.readlines()
    fin.close()
    remove_tmp_file('tmp.break')
    atm_breaks = []
    atm_start = 0
    for idx,line in enumerate(lines):
        line = line.strip().split()
        atm_break = [atm_start,int(line[0])]
        atm_breaks.append(atm_break)
        atm_start = int(line[0])

    return atm_breaks


def molden_gto_reader(natm,coords,fgto):
    
    fin = open(fgto,'r')
    lines = fin.readlines()
    fin.close()
    lines_iter = iter(lines)
    cent = []
    basis = {}

    def read_one_bas(lsym, ncont,fac=1):
        fac = float(fac)
        bas = []
        for i in range(int(ncont)):
            dat = (next(lines_iter)).split()
            bas.append([str(lsym).upper(),int(ncont),i,float(dat[0]), float(dat[1])*fac])
        return bas

    for line in lines_iter:
        dat = line.split()
        if len(dat) == 0:
            continue
        if dat[0].isdigit():
            iatm = int(dat[0]) - 1
            coord = coords[iatm]
            cent.append([iatm,coord])
            basis[iatm] = []
            #print("cent= ",cent)
        elif dat[0].upper() in 'SPDFGHIJ':
            #print("*dat= ", *dat)
            basis[iatm].append(read_one_bas(*dat))

    #print("cent[0]= ",cent[0])
    #print("len(basis[0])= ", len(basis[0]))
    #for i in range(len(basis[0])):
    #    print(i,basis[0][i])

    return cent,basis

def ebasis_writer(natm,centL,basis,fileout):
    fout = open(fileout,'w')
    buff = 'NO. OF UNIQUE CENTERS' + '\n'
    buff += str(natm) + '\n'
    fout.write(buff)
    for iatm in range(natm):
        icent = centL[iatm][1]
        #print("icent= ",icent)
        buff = "---------------"+icent[1]+"---------------" + '\n'
        icoord = icent[2]
        s = dp_format(icoord[0],icoord[1],icoord[2])
        buff += "CENTER" + '\n' + s + '\n'
        buff += 'SYMBOLS' + '\n' +"  "+ str(len(basis[iatm])) + '\n'
        fout.write(buff)
        for iprim in range(len(basis[iatm])):
            sym     = basis[iatm][iprim][0][0]
            ncont   = basis[iatm][iprim][0][1]
            #print("sym= ",sym)
            #print("ncont= ",ncont)
            buff  = sym + '\t' + str(ncont) + '\n'
            fout.write(buff)
            for icont in range(ncont):
                alpha = basis[iatm][iprim][icont][3]
                coeff = basis[iatm][iprim][icont][4]
                #print("alpha= ",alpha)
                #print("coeff= ",coeff)
                buff  = '\t'+ str(icont+1) + '\t' + dp_format(alpha) + '\t' + dp_format(coeff) + '\n'
                fout.write(buff)

    fout.close()
    return None


def read_coords(natm,fcoord):
    ang_to_bohr =1.8897259886e0
    fin = open(fcoord,'r')
    coords = []
    for i in range(natm):
        line = fin.readline()
        line = line.strip().split()
        atm = str(line[0])
        iatm = int(line[1])
        inuc = int(line[2])
        x    = float(line[3]) * ang_to_bohr
        y    = float(line[4]) * ang_to_bohr
        z    = float(line[5]) * ang_to_bohr
        xyz  = [x,y,z]
        coords.append([iatm,atm,xyz])
    fin.close()
    return coords

def read_coord_v2(natm,fcoord):
    ang_to_bohr =1.8897259886e0
    fin = open(fcoord,'r')
    coords = []
    xmin = np.zeros(3)
    xmax = np.zeros(3)
    for i in range(3):
        xmin[i] = 1.0e6
        xmax[i] = -1.0e6
    for i in range(natm):
        line = fin.readline()
        line = line.strip().split()
        atm = str(line[0])
        iatm = int(line[1])
        inuc = int(line[2])
        x    = float(line[3]) * ang_to_bohr
        y    = float(line[4]) * ang_to_bohr
        z    = float(line[5]) * ang_to_bohr
        xyz  = [x,y,z]
        for idim in range(3):
            if xmin[idim] > xyz[idim]:
                xmin[idim] = xyz[idim]
            if xmax[idim] < xyz[idim]:
                xmax[idim] = xyz[idim]
        coords.append([iatm,atm,xyz])
    fin.close()
    return xmin,xmax,coords


def test_momat(nb,roots,momat):
    fout = open('test_momat.inp','w')
    for i in range(nb):
        buff = "Ene= "+dp_format(roots[i]) +'\n' + 'SPIN' + '\n' + 'OCCUP' + '\n'
        fout.write(buff)
        for j in range(nb):
            buff = '\t'+ str(j+1)+ ' '+dp_format(momat[i,j]) + '\n'
            fout.write(buff)
    fout.close()
    return None


def parse(filename,outtag):
    print("----Parsing "+filename+"----",flush=True)
    O = TCM(outdir=outtag)

    extract_lines_from_file(O.coord_st, O.coord_ed,filename,O.fcoord)
    extract_lines_from_file(O.coord_ed,O.gto_end,filename,O.fgto)
    extract_lines_from_file(O.gto_end,O.mo_end,filename,O.fmovec,EOF_ext=False)
    grep_str_to_file(O.root_str,O.fmovec,O.froots)
    nb = line_cnt_file(O.froots)
    #roots,momat = extract_momat_molden(nb,nb,O.fmovec)
    spin,occ,roots,momat = extract_momat_v2(nb,nb,O.fmovec)
    write_mat(nb,momat,O.fmat)
    write_vec(nb,roots,O.froot)
    write_vec(nb,spin,O.fspin)
    write_vec(nb,occ,O.focc)

    natm = line_cnt_file(O.fcoord)
    xmin,xmax,coords = read_coord_v2(natm,O.fcoord)
    centL, basis = molden_gto_reader(natm,coords,O.fgto)
    ebasis_writer(natm,centL,basis,O.febasis)

    return xmin,xmax,nb

def write_grid_maker(xmin,xmax,nb,filename):
    fout = open(filename,'w')
    buff = 'NDIM' + '\n' + str(3) + '\n'
    buff += 'IDIM  NPT_IDIM IMIN      IMAX' + '\n'
    for i in range(3):  
        buff += str(i+1) + '\t' + str(100) + '\t' + dp_format(xmin[i]) + '\t' + dp_format(xmax[i]) + '\n'
    buff += 'NMO_for_grid NMO_MAX' + '\n'
    buff += 'XXX' + '\t' + str(nb) + '\n'
    buff += 'IDX IMO' + '\n'
    fout.write(buff)
    fout.close()
    return None

def test1A():
    #filename is the molden file path#
    #outtag is the directory path that you wish to write all outfiles to"
    filename = './water.molden'
    outtag   = './'
    xmin,xmax,nb = parse(filename,outtag)
    print("Finished parsing")
    print("xmin= ",xmin)
    print("xmax= ",xmax)
    delta = 12.0e0
    xmin = xmin - delta
    xmax = xmax + delta
    fgrid = outtag + 'grid_maker_info.inp'
    write_grid_maker(xmin,xmax,nb,fgrid)
    return None

if __name__ == "__main__" :
    assert sys.version_info[0] <= 3,"---MUST USE PYTHON 3.x---"
    now = datetime.datetime.now()
    print("START DATE: ",now)
    print()
    test1A()
    now = datetime.datetime.now()
    print("END DATE: ",now)
