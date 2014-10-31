import sys

def init(args):
    try:
        import __main__
        __main__.pymol_argv=['pymol',args]
        import pymol
        pymol.finish_launching()
    except:
        print "Unable to launch pymol"
        sys.exit()
    return pymol