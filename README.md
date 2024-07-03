# New_caledonia_RNAseq
```bash
# orthologous detection
# (base) kang1234@celia-PowerEdge-T640 Fri Jun 28 15:31:02 ~/New_caledonia/kraken
mkdir orthologue
cp Clun_runDrap/e-rmbt_editing/all_contigs.second_pass.fa orthologue/Clun_assembly.fa
cp Daru_runDrap/e-rmbt_editing/all_contigs.second_pass.fa orthologue/Daru_assembly.fa
cp Zlep_runDrap/e-rmbt_editing/all_contigs.second_pass.fa orthologue/Zlep_assembly.fa
# 然后用TransDecoder得到其蛋白质序列
nohup TransDecoder.LongOrfs -t Clun_assembly.fa  > Clun_transdecoder.process 2>&1 &
# [1] 27672
nohup TransDecoder.LongOrfs -t Daru_assembly.fa > Daru_transdecoder.process 2>&1 &
# [2] 27707
nohup TransDecoder.LongOrfs -t Zlep_assembly.fa > Zlep_transdecoder.process 2>&1 &
# [3] 27726

# (base) kang1234@celia-PowerEdge-T640 Fri Jun 28 16:04:37 ~/New_caledonia/kraken/orthologue
mkdir orthofinder_input_pep/
cd orthofinder_input_pep/
# (base) kang1234@celia-PowerEdge-T640 Fri Jun 28 16:07:12 ~/New_caledonia/kraken/orthologue/orthofinder_input_pep
less ../Clun_assembly.fa.transdecoder_dir/longest_orfs.pep|perl -alne 'if (/>/){$i++;my $nm=Clun."_"."$i";print">$nm"}else{print"$_"}' > Clun_pep.fasta
less ../Daru_assembly.fa.transdecoder_dir/longest_orfs.pep|perl -alne 'if (/>/){$i++;my $nm=Daru."_"."$i";print">$nm"}else{print"$_"}' > Daru_pep.fasta
less ../Zlep_assembly.fa.transdecoder_dir/longest_orfs.pep|perl -alne 'if (/>/){$i++;my $nm=Zlep."_"."$i";print">$nm"}else{print"$_"}' > Zlep_pep.fasta
cd ../
# (base) kang1234@celia-PowerEdge-T640 Fri Jun 28 16:07:13 ~/New_caledonia/kraken/orthologue/
nohup orthofinder -f orthofinder_input_pep -a 30 >orthofinder-process 2>&1 &
```
