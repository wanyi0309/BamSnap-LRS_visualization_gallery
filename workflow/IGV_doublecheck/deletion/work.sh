# At some loci, where BamSnap-LRS did not clearly display the supporting SV reads alignment, similar limitations were observed in SVhawkeye, Wally
# DEL24	chr5	626,820	629,092
# DEL40	chr11	99,320,935	99,323,166
# DEL41	chrY	51,102,780	51,118,029
# DEL42	chr6	104,464,281	104,490,692
# DEL43	chr1	227,886,102	227,898,290
# DEL44	chr4	67,265,086	67,286,141
# DEL45	chr1	227,872,682	227,889,341
# DEL47	chr4	71,956,601	72,068,816
# DEL48	chr6	26,614,353	26,637,304
# DEL49	chr6	26,604,794	26,629,112
# DEL50	chr5	21,488,029	21,568,569
# Therefore, we extracted the reads from these loci and took screenshots in IGV for verification. 
# For instructions on the BAM files, please refer to https://github.com/wanyi0309/Bamsnap-LRS_visualization_gallery/blob/main/workflow/gallay_used/deletion_work.sh

samtools view -h HIFI.primary.sort.bam chr5:626820-629092 chr11:99320935-99323166 chrY:51102780-51118029 chr6:104464281-104490692 chr1:227886102-227898290 chr4:67265086-67286141 chr1:227872682-227889341 chr4:71956601-72068816 chr6:26614353-26637304 chr6:26604794-26629112 chr5:21488029-21568569 | samtools sort -o DEL.HIFI.select.bam -
samtools index DEL.HIFI.select.bam

samtools view -h ONT.primary.sort.bam chr5:626820-629092 chr11:99320935-99323166 chrY:51102780-51118029 chr6:104464281-104490692 chr1:227886102-227898290 chr4:67265086-67286141 chr1:227872682-227889341 chr4:71956601-72068816 chr6:26614353-26637304 chr6:26604794-26629112 chr5:21488029-21568569 | samtools sort -o DEL.ONT.select.bam -
samtools index DEL.ONT.select.bam
