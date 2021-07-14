#_____ HELPER SCRIPTS _______________________________________________________#
def check_copy_finished(wc):
    return(["run_data/" + x + "/copy_finished" for x in map_samples_run[wc.sample]])

def get_fastqs(wc):
    folders = []
    if (wc.sample in map_samples_barcodes):
        folders.extend(["run_data/" + x  for x in map_samples_run[wc.sample]])
    elif (config['use_failed_reads']):
        folders.extend(["run_data/" + x + "/fastq_fail" for x in map_samples_run[wc.sample]])
    else: 
        folders.extend(["run_data/" + x + "/fastq_pass" for x in map_samples_run[wc.sample]])
    return{'folders': folders}

def get_summary_files(wc):
    g = "Sample_" + wc.sample + '/**/sequencing_summary*'
    files = [str(f) for f in glob(g, recursive=True)]
    #if (wc.sample in map_samples_barcodes):
    #    files.extend([
    #    "Sample_" +
    #    wc.sample + 
    #    "/runs/" + 
    #    r + 
    #    "sequencing_summary_split_" +
    #    wc_sample + 
    #    "_" +
    #    r +
    #    ".txt"
    #    for r in map_samples_run[wc.samples]])
    return{'summary_files': files}

def get_summary_file_split(wc):
    return

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

# Extract chromosome names from target region bed file
def get_chromosomes():
    with open(config['ref']['target_region']) as f:
        chrs = [row.split()[0] for row in f]
        return(chrs)

# Load additional project config files and overwrite default config
def load_project_config(f):
    if os.path.isfile(f):
        configfile: f
    else:
        print("Project config not available")
        pass
