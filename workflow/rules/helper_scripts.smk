#_____ HELPER SCRIPTS _______________________________________________________#
def get_input_folders(wc):
    folders = map_samples_folder[wc.sample]
    folders_exist = [x for x in folders if os.path.exists(x)]
    return{'folders': folders_exist} 

def get_summary_files(wc):
    g = "Sample_" + wc.sample + '/**/sequencing_summary*'
    files = [str(f) for f in glob(g, recursive=True)]
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
