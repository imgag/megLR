#_____ HELPER SCRIPTS _______________________________________________________#

def get_input_folders(wc):
    """
    Return the FASTQ data folder for a sample.
        Allows multiple runs per sample (restarted, repeated)
        Demultiplexed samples on a single flowcell with barcode information
        Includes folder with failed reads when specified in config
    """
    folders = map_samples_folder[wc.sample].copy()
    folders.append(['/fastq_pass/'.join(x) for x in map_samples_barcode[wc.sample]])
    print("Input Folders:" + str(folders[0]))
    if config['use_failed_reads']:
       folders.append(['/fastq_fail/'.join(x) for x in map_samples_barcode[wc.sample]])
    folders_exist = [x for x in folders[0] if os.path.exists(x)]
    return{'folders': folders[0]}

def lookup_split_summary_file(wc):
    """
    Looks up which split barcode file belongs to the sample in the wildcard.
    """
    bc=[unpack(x)[1] for x in map_samples_barcode[wc.sample]][0]
    #print(bc)
    f=[unpack(x)[0] for x in map_samples_barcode[wc.sample]][0]
    #print(f)
    split_folder = "qc/pycoqc/split_barcodes/sequencing_summary_"+bc+".txt" 
    return split_folder

def get_summary_files(wc):
    """
    Get all summary files that belong to a single sample.
    Requests summaries split by barcodes if necessary
    """
    folders = map_samples_folder[wc.sample].copy()
    files = [s for t in [glob(x+"**/sequencing_summary*", recursive = True) for x in folders] for s in t]
    if map_samples_barcode:
        folders_barcode = ['Sample_' + wc.sample for x in map_samples_barcode[wc.sample]]
        files += [x+"/sequencing_summary_bc_"+ wc.sample+".txt" for x in folders_barcode]
    
    return{'summary_files': files}

def aggregate_sample_pycoqc(wc):
    """
    Function that validates the checkpoint and checks for generated sample_pycoqcs
    that can be used in the multiqc report
    """
    barcode_qcs = []
    if ID_barcode_folders:
        [checkpoints.split_summary_perbarcode.get(folder=x).output[0] for x in ID_samples]
    else:
        return(barcode_qcs)

    #if map_samples_barcode: 
    #    for k,v in map_samples_barcode.copy().items():
    #        if k in ID_folders:
    #            print("Found" + k)
    #        folders_barcode = v[0]
    #        g = glob(os.path.join(checkpoint_output, "/summary_statistics_{bc}.txt"))
    #        qc_out_selected += (expand("qc/pycoqc/split_{folder}/summary_statistics_{bc}.txt",
    #            folder = folders_barcode,
    #            bc = g))


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
