"""
Here are just comments for now, might be implemented later stages

Using Guppy v4.0.11
Model:

'dna_r9.4.1_450bps_hac_prom'

guppy_basecaller \
    --input_path 20190425_I19D001a04Mega15_BP10 \
    --recursive \
    --save_path a04_01 \
    --device cuda:0 \
    --config dna_r9.4.1_450bps_hac_prom.cfg &

guppy_basecaller \
    --input_path 20190502_I19D001a04_Mega15_BP10_AMPure_1X \
    --recursive \
    --save_path a04_02 \
    --device cuda:1 \
    --config dna_r9.4.1_450bps_hac_prom.cfg &
"""