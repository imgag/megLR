#_____ HELPER SCRIPTS _______________________________________________________#
def check_copy_finished(wc):
    return(["run_data/" + x + "/copy_finished" for x in map_samples_run[wc.sample]])

def get_fastqs(wc):
    folders = []
    folders.extend(["run_data/" + x + "/fastq_pass" for x in map_samples_run[wc.sample]])
    if (config['use_failed_reads']):
        folders.extend(["run_data/" + x + "/fastq_fail" for x in map_samples_run[wc.sample]])
    return{'folders': folders}

def get_summary_files(wc):
    g = "Sample_" + wc.sample + '/**/sequencing_summary*'
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
    