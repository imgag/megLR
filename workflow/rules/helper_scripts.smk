#_____ HELPER SCRIPTS _______________________________________________________#

def get_input_folders(wc):
    """
    Return the FASTQ data folder for a sample.
        Allows multiple runs per sample (restarted, repeated)
        Demultiplexed samples on a single flowcell with barcode information
        Includes folder with failed reads when specified in config
    """
    folders = map_samples_folder[wc.sample].copy()
    folders += ['/pass/'.join(x) for x in map_samples_barcode[wc.sample]]
    if config['use_failed_reads']:
       folders += ['/fail/'.join(x) for x in map_samples_barcode[wc.sample]]
    folders_exist = [x for x in folders if os.path.exists(x)]
    return{'folders': folders_exist} 

def lookup_split_summary_file(wc):
    """
    Looks up which split barcode file belongs to the sample in the wildcard.
    """
    #files = [glob(x+"/sequencing_summary*.txt") for x in wc.folders]

    outfolder = checkpoints.split_summary_perbarcode.get(folder=wildcards.)
    f = [x[0] for x in map_samples_barcode[wc.sample]]
    bc = [x[0] for x in map_samples_barcode[wc.sample]]
    split_folders = ["qc/pycoqc/split_" + x for x in folders_barcode] 
    return "qc/pycoqc/split_21033_rebasecalled/"

def get_summary_files(wc):
    """
    Get all summary files that belong to a single sample.
    Requests summaries split by barcodes if necessary
    """
    folders = map_samples_folder[wc.sample].copy()
    files = [glob(x+"/sequencing_summary*.txt") for x in folders]
    folders_barcode = ['Sample_' + wc.sample for x in map_samples_barcode[wc.sample]]
    files += [x+"/sequencing_summary_bc_"+ wc.sample+".txt" for x in folders_barcode]
    print(files)
    #return "qc/pycoqc/split_21033_rebasecalled/"
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
