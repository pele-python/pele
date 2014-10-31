import os

def read_text_file(filename):
    '''
    Returns data read from a text file as an array of strings. 
    Returns false on failure.
    '''
    try:
        vst = open(filename, 'r')
        raw_data=vst.readlines()
        vst.close()
    except IOError:
        raw_data=False
        
    return raw_data

def write_text_file(filename, lines):
    '''
    Writes data read to a text file as an array of strings. 
    Returns false on failure.
    '''
    output=True
    try:
        vst = open(filename, 'w')
        vst.writelines(lines)
        vst.close()
    except IOError:
        output=False
      
    return output

def read_atoms_from_pdb(pdb_strings):
    '''
    Loops through an array of strings  - one for line in a pdb file. 
    Each ATOM or HETATM line is converted into an array of data.
    '''
    atoms = []
    for line in pdb_strings:
        if line[0:4]=="ATOM" or line[0:6]=="HETATM":
            try:
                a = parsePdbLine(line)
                if not a in atoms:
                    atoms.append(a)
            except:
                print "line: " + line + " not understood"
                exit(0)
    return atoms

def parsePdbLine(line):
    '''
    Extracts each piece of information a single ATOM or HETATM line
    from a PDB file.
    '''
    l=line
    atom_seri = int(l[6:11])
    atom_name = l[12:16].split()[0]
    alte_loca = l[16]
    resi_name = l[17:20].split()[0]
    chai_iden = l[21]
    resi_numb = int(l[22:26])
    code_inse = l[26]
    atom_xcoo = float(l[30:38])
    atom_ycoo = float(l[38:46])
    atom_zcoo = float(l[46:54])
    try:
        atom_occu = float(l[54:60])
    except:
        atom_occu=0.0

    try:
        atom_bfac = float(l[60:66])
    except:
        atom_bfac=0.0    
    
    try:
        seg_id = l[72:76]
    except:
        seg_id=' '

    try:
        atom_symb = l[76:78].split()[0]
    except:
        try:
            atom_symb = l[68]
        except:
            atom_symb= ' '

    try:
        charge=l[78:80]
    except:
        charge=' '

    return [atom_seri, atom_name, alte_loca, resi_name, chai_iden, resi_numb, code_inse, atom_xcoo, atom_ycoo, atom_zcoo, atom_occu, atom_bfac,seg_id,atom_symb,charge]
