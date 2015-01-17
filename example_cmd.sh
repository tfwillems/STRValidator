bams=`ls /Users/tfwillems/Desktop/Club_files/chry/final_callset/chrY_BAMs/filtered_bams/*.chrY.bam -m | tr -d '\n'`
bais=`ls /Users/tfwillems/Desktop/Club_files/chry/final_callset/chrY_BAMs/filtered_bams/*.chrY.bam.bai -m | tr -d '\n'`
fasta_dir="/Users/tfwillems/Desktop/STR_files/chromFa/"
vizalign="/Users/tfwillems/Desktop/Coding/HipSTR/vizalign"

./str_validator.py --bams $bams --bais $bais --fasta $fasta_dir --vizalign $vizalign
