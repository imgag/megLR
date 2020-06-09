#=============================================================================#
#                    R E A D    P R E P A R A T I O N                         #
#=============================================================================#


#_____ READ COPY FROM PROMETHION _____________________________________________#

rule copy_prom:
    input:
        config['sample_file']
    output:
        touch("run_data/{run}/copy_finished")
    log:
        "logs/{run}_copy.log"
    threads: 1
    shell:
        """
        folders=$(ssh microbio "find /mnt/promdata -type d -name \\*{wildcards.run}\\* 2> /dev/null; return 0;") > {log} 2>> {log}
        echo "Found the following folders: " >> {log} 2>> {log}
        echo "$folders" >> {log} 2>> {log}

        for folder in $folders
        do
            rsync \
                --archive \
                --inplace \
                --exclude="*.fast5" \
                --log-file {log} \
                --verbose \
                microbio:$folder $PWD/run_data/
        done
        """

#_____ COMPRESS TO SINGLE .FASTQ.GZ __________________________________________#

rule join_fastq:
    input:
        check_copy_finished
    output:
        "Sample_{sample}/{sample}.fastq.gz"
    threads:
        8
    shell:
        """
        find {input} -name '*.fastq' -exec cat {{}} + | pigz -p {threads} -c > {output}
        """