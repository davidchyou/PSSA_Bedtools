SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
INFILE=$1
MAX=100
$SCRIPTPATH/bedtools closest -a $INFILE -b $INFILE -s -iu -io -D a 2>/dev/null | \
awk -F'\t' '($13>0&&$13<='$MAX')' | \
perl -F'\t' -lane 'BEGIN{%chemical_names=();%init_seqs=()}
{($drug_A,$drug_A_name,$drug_B,$drug_B_name)=@F[(4,3,10,9)];
$chemical_names{$drug_A}=$drug_A_name;
$chemical_names{$drug_B}=$drug_B_name;
if (exists $init_seqs{$drug_A}{$drug_B})
{$init_seqs{$drug_A}{$drug_B}++;} 
else {$init_seqs{$drug_A}{$drug_B}=1;}
}END{@lst_drugs=keys(%chemical_names);
for($i=0;$i<scalar(@lst_drugs);$i++){
for($j=0;$j<scalar(@lst_drugs);$j++){
$a=0;$b=0;$ap=1;$bp=1;$aid=@lst_drugs[$i];$bid=@lst_drugs[$j];
$a_name=$chemical_names{$aid};$b_name=$chemical_names{$bid};
$a=(exists $init_seqs{$aid}{$bid})?($init_seqs{$aid}{$bid}):0;
$b=(exists $init_seqs{$bid}{$aid})?($init_seqs{$bid}{$aid}):0;
$ap=$a+1;$bp=$b+1;$log_sr=log($ap/$bp);$se=sqrt((1/$ap)+(1/$bp));
$ci_low=$log_sr-(1.96*$se);$ci_high=$log_sr+(1.96*$se);$signal=($ci_low>0)?1:0;
@rs=($aid,$bid,$a_name,$b_name,$a,$b,$ap,$bp,$log_sr,$ci_low,$ci_high,$signal);
$strout=join("\t",@rs);print $strout;
}}}'
