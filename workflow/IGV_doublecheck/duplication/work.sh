# At some loci, where BamSnap-LRS did not clearly display the supporting SV reads alignment, similar limitations were observed in SVhawkeye, Wally
# DUP1	chr20	62,753,581	62,757,947
# DUP2	chr5	14,632,470	14,636,733
# DUP3	chr4	83,000,459	83,004,511
# DUP4	chr4	113,886,291	113,890,691
# DUP5	chr5	1,096,166	1,100,439
# DUP6	chr5	70,985	75,527
# DUP7	chr7	150,033,595	150,037,935
# DUP10	chr12	131,328,035	131,332,271
# DUP11	chr22	45,226,893	45,231,369
# DUP12	chr8	99,032,443	99,036,973
# DUP15	chr14	105,626,058	105,630,340
# Therefore, we extracted the reads from these loci and took screenshots in IGV for verification. 
# For instructions on the BAM files, please refer to https://github.com/wanyi0309/Bamsnap-LRS_visualization_gallery/blob/main/workflow/gallay_used/insertion_work.sh

samtools view -h -M SRR16005301.ont.hg38.sorted.bam chr20:62753581-62757947 chr5:14632470-14636733 chr4:83000459-83004511 chr4:113886291-113890691 chr5:1096166-1100439 chr5:70985-75527 chr7:150033595-150037935 chr12:131328035-131332271 chr22:45226893-45231369 chr8:99032443-99036973 chr14:105626058-105630340 | samtools sort -o DUP.hg38.ONT.select.bam -
samtools index DUP.hg38.ONT.select.bam

samtools view -h -M SRR8955953.pacbio.hg38.sorted.bam chr20:62753581-62757947 chr5:14632470-14636733 chr4:83000459-83004511 chr4:113886291-113890691 chr5:1096166-1100439 chr5:70985-75527 chr7:150033595-150037935 chr12:131328035-131332271 chr22:45226893-45231369 chr8:99032443-99036973 chr14:105626058-105630340 | samtools sort -o DUP.hg38.Pacbio.select.bam -
samtools index DUP.hg38.Pacbio.select.bam
