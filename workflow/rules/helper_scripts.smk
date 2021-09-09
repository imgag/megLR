#_____ HELPER SCRIPTS _______________________________________________________#

def get_input_folders(wc):
    """
    Return the FASTQ data folder for a sample.
        Allows multiple runs per sample (restarted, repeated)
        Demultiplexed samples on a single flowcell with barcode information
        Includes folder with failed reads when specified in config
    """
    folders = ['/pass/'.join(x) if (type(x) == tuple) else x for x in  map_samples_folder[wc.sample]]
    if config['use_failed_reads']:
       folders.append(['/fail/'.join(x) for x in map_samples_folder[wc.sample] if (type(x) == tuple)])
    folders_exist = [x for x in folders if os.path.exists(x)]
    return{'folders': folders_exist} 

def get_summary_files(wc):
    """
    Get the summary files for a single sample.
    Let PycoQC Split the aggregated summary file by barcodes if necessary
    """
    #TODO UNFINISHED WORK
    folders = map_samples_folder[wc.sample]
    barcoded = [x for x in folders if type(x) == tuple]
    files = [str(f) for f in glob(g, recursive=True)]
    return{'summary_files': files}

def print_message():
    print('')
    print('  ██████╗ ███╗   ██╗████████╗  ████████╗ ██████╗  ██████╗ ██╗     ███████╗')
    print(' ██╔═══██╗████╗  ██║╚══██╔══╝  ╚══██╔══╝██╔═══██╗██╔═══██╗██║     ██╔════╝')
    print(' ██║   ██║██╔██╗ ██║   ██║        ██║   ██║   ██║██║   ██║██║     ███████╗')
    print(' ██║   ██║██║╚██╗██║   ██║        ██║   ██║   ██║██║   ██║██║     ╚════██║')
    print(' ╚██████╔╝██║ ╚████║   ██║███████╗██║   ╚██████╔╝╚██████╔╝███████╗███████║')
    print('  ╚═════╝ ╚═╝  ╚═══╝   ╚═╝╚══════╝╚═╝    ╚═════╝  ╚═════╝ ╚══════╝╚══════╝') 
    print()
    print('Pipeline runs these steps:')
    print(*config['steps'], sep=" | ")
    pass

def get_chromosomes():
    """Extract chromosome names from target region bed file"""
    with open(config['ref']['target_region']) as f:
        chrs = [row.split()[0] for row in f]
        return(chrs)

def load_project_config(f):
    """Load additional project config files and overwrite default config"""
    if os.path.isfile(f):
        configfile: f
    else:
        print("Project config not available")
        pass
