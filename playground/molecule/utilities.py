import os

def read_text_file(filename):
    '''
    Returns data read from a text file as an array of strings. 
    Returns false on failure.
    '''
    try:
        f = open(filename, 'r')
        raw_data=f.readlines()
        f.close()
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
        f = open(filename, 'w')
        f.writelines(lines)
        f.close()
    except IOError:
        output=False
      
    return output