#_____ HELPER SCRIPTS _______________________________________________________#
def check_copy_finished(wc):
    return(["run_data/" + x + "/copy_finished" for x in map_samples_run[wc.sample]])

def get_summary_files(wc):
    g = "Sample_" + wc.sample + '/**/sequencing_summary*'
    files = [str(f) for f in glob(g, recursive=True)]
    return{'summary_files': files}
