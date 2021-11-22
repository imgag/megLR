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
    bc=[unpack(x)[1] for x in map_samples_barcode[wc.sample]][0]
    #print(bc)
    f=[unpack(x)[0] for x in map_samples_barcode[wc.sample]][0]
    #print(f)
    split_folder = "qc/pycoqc/split_" + f +"/sequencing_summary_"+bc+".txt" 
    return split_folder

def get_summary_files(wc):
    """
    Get all summary files that belong to a single sample.
    Requests summaries split by barcodes if necessary
    """
    folders = map_samples_folder[wc.sample].copy()
    files = [s for t in [glob(x+"/**/sequencing_summary*") for x in folders] for s in t]
    if map_samples_barcode:
        folders_barcode = ['Sample_' + wc.sample for x in map_samples_barcode[wc.sample]]
        files += [x+"/sequencing_summary_bc_"+ wc.sample+".txt" for x in folders_barcode]
    
    return{'summary_files': files}

def aggregate_multiqc_input(wc):
    """
    Function that validates the checkpoint and checks for generated sample_pycoqcs
    that can be used in the multiqc report
    """

    qc_out = {
    'mapping' : expand("qc/qualimap/{s}_genome/genome_results.txt", s = ID_samples),
    'assembly' : ["qc/quast_results/report.tsv"],
    'variant_calling':[expand("qc/variants/{s}.stats", s = ID_samples)],
    'structural_variant_calling' : [],
    'cDNA_stringtie' : expand("qc/gffcompare/{s}_stringtie/{s}_stringtie.stats", s = ID_samples) +
        expand("qc/pychopper/{s}_stats.txt", s = ID_samples), 
    'cDNA_flair': 
        expand("qc/rseqc/{s}.read_distribution.txt", s = ID_samples) + 
        expand("qc/rseqc/{s}.geneBodyCoverage.txt", s = ID_samples) +    
        expand("qc/gffcompare/{s}_flair/{s}_flair.stats", s = ID_samples),
    'cDNA_expression' : 
        #expand("qc/qualimap/{s}_rna/rnaseq_qc_results.txt", s = ID_samples) + 
        expand("qc/rseqc/{s}.read_distribution.txt", s = ID_samples) + 
        expand("qc/rseqc/{s}.geneBodyCoverage.txt", s = ID_samples) + 
        expand("qc/pychopper/{s}_stats.txt", s = ID_samples) +
        expand("Sample_{s}/{s}.counts.tsv.summary", s = ID_samples),
    'cDNA_pinfish' : [],
    'dual_demux' : [],
    'de_analysis' : [],
    'qc' : ["qc/pycoqc/per_run/run_multiqc_report.html",
        expand("qc/pycoqc/per_sample/{s}.pycoQC.json", s = ID_samples)],
    }

    # Additional output options
    if config['vc']['create_benchmark']:
        qc_out['variant_calling'] +=  expand("qc/happy/{s}.summary.csv", s=ID_samples)

    qc_out_selected = [qc_out[step] for step in config['steps']]    

    if map_samples_barcode: 
        print(map_samples_barcode)
        print(len(map_samples_barcode.items()))
        for k,v in map_samples_barcode.copy().items():
            print(k)
            print(v)
            folders_barcode = v[0][0]
            print(folders_barcode)
            checkpoint_output=checkpoints.split_summary_perbarcode.get(folder=folders_barcode[0]).output[0]
            g = glob(os.path.join(checkpoint_output, "/summary_statistics_{bc}.txt"))
            qc_out_selected += (expand("qc/pycoqc/split_{folder}/summary_statistics_{bc}.txt",
                folder = folders_barcode,
                bc = g))
    
    return(qc_out_selected)


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
