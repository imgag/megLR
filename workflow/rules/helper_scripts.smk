#_____ HELPER SCRIPTS _______________________________________________________#

def get_input_folders(wc):
    """
    Return the FASTQ data folder for a sample.
        Allows multiple runs per sample (restarted, repeated)
        Demultiplexed samples on a single flowcell with barcode information
        Includes folder with failed reads when specified in config

        New option: Ifq folder contains a subfolder called "fastq_rebasecalled", 
        reads will be taken only from this folder. Option can be turned off in config option
        'fastq_prefer_rebasecalled' 
    """
    folders = map_samples_folder[wc.sample].copy()
    if config['verbose']: print("Input Folders:" + str(folders), end = '')
    
    if map_samples_barcode:
        [folders.append('/fastq_pass/'.join(x)) for x in map_samples_barcode[wc.sample]]
        if config['use_failed_reads']:
            folders.append(['/fastq_fail/'.join(x) for x in map_samples_barcode[wc.sample]])
    
    folders_updated = list()

    if config['fastq_prefer_rebasecalled']:
        for f in folders:
            f_rebasecalled =  glob(f+"/**/fastq_rebasecalled", recursive = True)
            if f_rebasecalled:
                [folders_updated.append(f) for f in f_rebasecalled]
            else:
                folders_updated.append(f)
    
    if config['verbose']: print(" | Updated to:" + str(folders_updated))
    return{'folders': folders_updated}

def get_input_folders_fast5(wc):
    """
    Return the FAST5 data folder for a sample.
        Allows multiple runs per sample (restarted, repeated)
        Demultiplexed samples on a single flowcell with barcode information
        Includes folder with failed reads when specified in config
    """
    folders = map_samples_folder[wc.sample].copy()
    folders.append(['/fast5_pass/'.join(x) for x in map_samples_barcode[wc.sample]])
    if config['verbose']: print("Input Folders:" + str(folders[0]))
    if config['use_failed_reads']:
       folders.append(['/fast5_fail/'.join(x) for x in map_samples_barcode[wc.sample]])
    folders_exist = [x for x in folders[0] if os.path.exists(x)]
    return{'folders': folders[0]}

def use_bam(wc):
    """
    Decide whether to use normal (Guppy) aligned bam
    or Bonito bam (with modified bases)
    """
    bam = "Sample_{s}/{s}.mod.bam".format(s=wc.sample) if config['use_mod_bam'] else "Sample_{s}/{s}.bam".format(s=wc.sample)
    return(bam)

def lookup_split_summary_file(wc):
    """
    Looks up which split barcode file belongs to the sample in the wildcard.
    """
    
    bc=[unpack(x)[1] for x in map_samples_barcode[wc.sample]][0]
    print("BC: " + bc)
    folder=[unpack(x)[0] for x in map_samples_barcode[wc.sample]][0]
    print("Folder: " + folder)
    for r, f in map_runs_folder.items():   
        if f and f[0] == folder:
            run1 = r
    print("Map Runs:", map_runs_folder)
    split_output_folder = checkpoints.split_summary_perbarcode.get(run = run1).output
    split_summary_file = os.path.join(split_output_folder, "sequencing_summary_"+bc+".txt")
    return split_summary_file

def get_summary_files_sample(wc):
    """
    Get all summary files that belong to a single sample.
    Requests summaries split by barcodes if necessary
    """
   
    if map_samples_barcode:
        files = ["Sample_{s}/sequencing_summary_bc_{s}.txt".format(s = wc.sample)]
   
    else:
        folders = map_samples_folder[wc.sample]
        files = [s for t in [glob(x+"/**/sequencing_summary*", recursive = True) for x in folders] for s in t]
        run=[r for f, r in map_runs_folder.items() if folder in f]

        if config['fastq_prefer_rebasecalled']:
            folders_rebasecalled = [s for t in [glob(x+"/**/fastq_rebasecalled", recursive = True) for x in folders] for s in t]
            if folders_rebasecalled:
                folders = folders_rebasecalled
                files = [s for t in [glob(x+"/**/sequencing_summary*", recursive = True) for x in folders_rebasecalled] for s in t]
    
    return {'summary_files': files}


def get_summary_files_to_split(wc):
    """
    Returns the summary statistic files to be split. Should return multiple files if same barcoded library was sequenced
    on multiple flowcells
    """

    files = []
    folders = []
    for r, f in map_runs_folder.items():
        if r == wc.run:
            folders.append(f[0])
    files = [s for t in [glob(x+"/**/sequencing_summary*", recursive = True) for x in folders] for s in t]
    return files


def get_db_report(wc):
    """
    Get Report .pdf and Report.md files
    Replace with empty dummy files if not available
    """
    reports = [x for y in [glob(r + "/**/report_*") for r in map_runs_folder[wc.run]] for x in y]
    md = [x for x in reports if x.endswith('.md')]
    reports = [x for x in reports if x not in md]

    #print(md)
    #print(reports)

    if not md:
        if config['verbose']: print("Warning: No markdown report found for run(s) " + wc.run)
        md = str(os.path.join(workflow.basedir, "../resources/dummyfiles/report.md"))
    
    return{'md': md, 'reports': reports}
    
def get_db_mux(wc):
    """
    Get mux stats file
    Replace with empty dummy files if not available
    """

    mux = [x for y in [glob(r + "/**/mux_scan_data*.csv") for r in map_runs_folder[wc.run]] for x in y]
    mux += [x for y in [glob(r + "/**/other_reports/pore_scan_data*.csv") for r in map_runs_folder[wc.run]] for x in y]
    
    if not mux:
        if config['verbose']: print("Warning: No mux stats (.csv) found for run(s) " + wc.run)
        mux = str(os.path.join(workflow.basedir, "../resources/dummyfiles/mux.csv"))
    
    return(mux)

def get_db_barcode(wc):
    """
    Get barcode tsv file
    Replace with empty dummy files if not available
    """

    bc = [x for y in [glob(r + "/**/barcode_alignment*.tsv") for r in map_runs_folder[wc.run]] for x in y]

    if not bc:
        if config['verbose']: print("Warning: No barcode file (.csv) found for run(s) " + wc.run)
        bc = str(os.path.join(workflow.basedir, "../resources/dummyfiles/barcodes.tsv"))
    
    return(bc)



def get_cnvkit_bam(wc):
    """
    Get bamfile for cnvkit either from config file (reference) or project folder (sample)
    """
    if wc.sample in config['cnvkit']['reference_samples']:
        bam = config['cnvkit']['reference_samples'][wc.sample]
    else:
        bam = "Sample_{s}/{s}.bam".format(s = wc.sample)
    if config['verbose']: print(bam)
    return{'bam': bam}

def get_cdna_bam(wc):
    """
    Return deduplicated or non-deduplicated cDNA BAM, depending on config options
    """
    if config['cdna']['with_umi']:
        return("Sample_{}/{}.spliced.dedup.bam".format(wc.sample))
    else:
        return("Sample_{}/{}.spliced.bam".format(wc.sample))

def print_message():
    print("                 _    ___ ")
    print("  _ __  ___ __ _| |  | _ \\")
    print(" | '  \/ -_) _` | |__|   /")
    print(" |_|_|_\___\__, |____|_|_\\")
    print("           |___/           ")
    print('')
    print('Pipeline branches:')
    print('╦')
    print(*["╠═ "+s for s in config['steps']], sep="\n")
    print('')
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
