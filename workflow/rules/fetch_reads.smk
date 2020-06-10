#_____ COPY READS FROM FROM PROMETHION _______________________________________#

rule copy_prom:
    input:
        config['sample_file']
    output:
        touch("run_data/{run}/copy_finished")
    log:
        "logs/{run}_copy.log"
    threads: 1
    params:
        exclude = "" if config['transfer_fast5'] else "*.fast5",
        ssh_profile = config['prom']['ssh_profile'],
        data_folder = config['prom']['data_folder']
    shell:
        """
        folders=$(ssh {params.ssh_profile} \
            "find {params.data_folder} -type d -name \\*{wildcards.run}\\* 2> /dev/null; return 0;" \
            ) > {log} 2>> {log}
        echo "Found the following folders: " >> {log} 2>&1
        echo "$folders" >> {log} 2>&1

        for folder in $folders
        do
            rsync \
                --archive \
                --inplace \
                --exclude="{params.exclude}" \
                --log-file {log} \
                --verbose \
                microbio:$folder $PWD/run_data/ \
                >> {log} 2>&1
        done
        """

#_____ COMPRESS TO SINGLE .FASTQ.GZ __________________________________________#

rule join_fastq:
    input:
        check_copy_finished,
        unpack(get_fastqs)
    output:
        "Sample_{sample}/{sample}.fastq.gz"
    threads:
        8
    conda:
        "../env/pigz.yml"
    shell:
        """
        find {input.folders} -name '*.fastq' -exec cat {{}} + \
            | pigz -p {threads} -c > {output}
        """
