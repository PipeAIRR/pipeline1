$HOSTNAME = ""
params.outdir = 'results'  

// Add for each process an option to change the parameters. Default is the set params
//* params.nproc =  1  //* @input @description:"How many processes to use for each step. Default 1"
//* params.edit_filter_quality_params =  "no"  //* @dropdown @options:"yes","no"  @show_settings:"filter_seq"
//* params.edit_mask_primer_params =  "no"  //* @dropdown @options:"yes","no" @show_settings:"MaskPrimers"
//* params.edit_pair_sequence_pre_consensus_params =  "no"  //* @dropdown @options:"yes","no" @show_settings:"pair_seq" @description:"If edit true, then parametrs in the pre_consensus tab should be edited"
//* params.edit_build_consensus_params =  "no"  //* @dropdown @options:"yes","no" @show_settings:"build_consensus"
//* params.edit_pair_sequence_post_consensus_params =  "no"  //* @dropdown @options:"yes","no" @show_settings:"pair_seq" @description:"If edit true, then parametrs in the post_consensus tab should be edited"
//* params.edit_assemble_pairs_params =  "no"  //* @dropdown @options:"yes","no" @show_settings:"assemble_pairs"
//* params.edit_parse_headers_params =  "no"   //* @dropdown @options:"yes","no" @show_settings:"parse_header" @description:"If edit true, then parametrs in the parse_header tab should be edited"
//* params.edit_collapse_seq_params =  "no"   //* @dropdown @options:"yes","no" @show_settings:"collapse_seq"
//* params.edit_split_seq_params =  "no"   //* @dropdown @options:"yes","no" @show_settings:"split_seq"
//* params.edit_parse_header_table =  "no"   //* @dropdown @options:"yes","no" @show_settings:"parse_headers" @description:"If edit true, then parametrs in the parse_header_table tab should be edited"

//* autofill
if ($HOSTNAME == "default"){
    $DOCKER_IMAGE = "immcantation/suite:4.3.0"
    $DOCKER_OPTIONS = "-v /work:/work"

}

//* platform
if ($HOSTNAME == "ig03.lnx.biu.ac.il"){
    $DOCKER_IMAGE = "immcantation/suite:4.3.0"
    $DOCKER_OPTIONS = "-v /work:/work"
	$CPU  = 48
    $MEMORY = 300 
}
//* platform


//* autofill

if((params.edit_filter_quality_params && (params.edit_filter_quality_params == "no"))){

    // Process Parameters for params.Filter_Sequence_Quality_filter_seq_quality:
    params.Filter_Sequence_Quality_filter_seq_quality.method = "quality"
    params.Filter_Sequence_Quality_filter_seq_quality.nproc = params.nproc
    params.Filter_Sequence_Quality_filter_seq_quality.q = "20"

}

if((params.edit_mask_primer_params && (params.edit_mask_primer_params == "no"))){
    // Process Parameters for Mask_Primer_MaskPrimers:
    params.Mask_Primer_MaskPrimers.method = ["score","score"]
    params.Mask_Primer_MaskPrimers.mode = ["cut","mask"]
    params.Mask_Primer_MaskPrimers.primer_field = ["PRIMER","PRIMER"]
    params.Mask_Primer_MaskPrimers.barcode_field = ["BARCODE","BARCODE"]
    params.Mask_Primer_MaskPrimers.start = [0,0]
    params.Mask_Primer_MaskPrimers.barcode = ["true","false"]
    params.Mask_Primer_MaskPrimers.umi_length = [15,0]
    params.Mask_Primer_MaskPrimers.maxerror = [0.2,0.2]
    params.Mask_Primer_MaskPrimers.revpr = ["false","false"]
    params.Mask_Primer_MaskPrimers.failed = "true"
}


if((params.edit_pair_sequence_pre_consensus_params && (params.edit_pair_sequence_pre_consensus_params == "no"))){

    // Process Parameters for params.Pair_Sequence_pre_consensus_pair_seq:
    params.Pair_Sequence_pre_consensus_pair_seq.coord = "sra"
    params.Pair_Sequence_pre_consensus_pair_seq.act = "none"
    params.Pair_Sequence_pre_consensus_pair_seq.copy_fields_1 = "BARCODE"
    params.Pair_Sequence_pre_consensus_pair_seq.copy_fields_2 = ""
    params.Pair_Sequence_pre_consensus_pair_seq.nproc = params.nproc

}


if((params.edit_build_consensus_params && (params.edit_build_consensus_params == "no"))){

    // Process Parameters for params.Build_Consensus_build_consensus:
    params.Build_Consensus_build_consensus.failed = "false"
    params.Build_Consensus_build_consensus.nproc = params.nproc
    params.Build_Consensus_build_consensus.barcode_field = ["BARCODE","BARCODE"]
    params.Build_Consensus_build_consensus.primer_field = ["PRIMER","PRIMER"]
    params.Build_Consensus_build_consensus.act = ["none","none"]
    params.Build_Consensus_build_consensus.mincount = [1,1]
    params.Build_Consensus_build_consensus.minqual = [0,0]
    params.Build_Consensus_build_consensus.minfreq = [0.6,0.6]
    params.Build_Consensus_build_consensus.maxerror = [0.1,0.1]
    params.Build_Consensus_build_consensus.prcons = [0.6,"none"]
    params.Build_Consensus_build_consensus.maxgap = [0.5,0.5]
    params.Build_Consensus_build_consensus.maxdiv = ["none","none"]
    params.Build_Consensus_build_consensus.dep = ["false","false"]

}

if((params.edit_pair_sequence_post_consensus_params && (params.edit_pair_sequence_post_consensus_params == "no"))){

    // Process Parameters for params.Pair_Sequence_post_consensus_pair_seq:
    params.Pair_Sequence_post_consensus_pair_seq.coord = "presto"
    params.Pair_Sequence_post_consensus_pair_seq.act = "none"
    params.Pair_Sequence_post_consensus_pair_seq.copy_fields_1 = ""
    params.Pair_Sequence_post_consensus_pair_seq.copy_fields_2 = ""
    params.Pair_Sequence_post_consensus_pair_seq.nproc = params.nproc

}

if((params.edit_assemble_pairs_params && (params.edit_assemble_pairs_params == "no"))){

    // Process Parameters for params.Assemble_pairs_assemble_pairs:
    params.Assemble_pairs_assemble_pairs.method = "align"
    params.Assemble_pairs_assemble_pairs.coord = "presto"
    params.Assemble_pairs_assemble_pairs.rc = "tail"
    params.Assemble_pairs_assemble_pairs.head_fields_R1 = "CONSCOUNT"
    params.Assemble_pairs_assemble_pairs.head_fields_R2 = "CONSCOUNT PRCONS"
    params.Assemble_pairs_assemble_pairs.failed = "false"
    params.Assemble_pairs_assemble_pairs.fasta = "false"
    params.Assemble_pairs_assemble_pairs.nproc = params.nproc
    params.Assemble_pairs_assemble_pairs.alpha = 0.00001
    params.Assemble_pairs_assemble_pairs.maxerror = 0.3
    params.Assemble_pairs_assemble_pairs.minlen = 8
    params.Assemble_pairs_assemble_pairs.maxlen = 1000
    params.Assemble_pairs_assemble_pairs.scanrev = "false"
    params.Assemble_pairs_assemble_pairs.minident = 0.5
    params.Assemble_pairs_assemble_pairs.evalue = 0.00001
    params.Assemble_pairs_assemble_pairs.maxhits = 100
    params.Assemble_pairs_assemble_pairs.fill = "false"
    params.Assemble_pairs_assemble_pairs.gap = 0
}


if((params.edit_parse_headers_params && (params.edit_parse_headers_params == "no"))){
    // Process Parameters for params.Parse_header_parse_headers:
    params.Parse_header_parse_headers.method = "collapse"
    params.Parse_header_parse_headers.act = "min"
    params.Parse_header_parse_headers.args = "-f CONSCOUNT"
}

if((params.edit_collapse_seq_params && (params.edit_collapse_seq_params == "no"))){
    // Process Parameters for params.edit_collapse_seq_params:
    params.collapse_sequences_collapse_seq.act = "sum"
    params.collapse_sequences_collapse_seq.max_missing = 20
    params.collapse_sequences_collapse_seq.inner = "true"
    params.collapse_sequences_collapse_seq.uf = "CREGION"
    params.collapse_sequences_collapse_seq.cf = "CONSCOUNT"
    params.collapse_sequences_collapse_seq.nproc = params.nproc
}

if((params.edit_split_seq_params && (params.edit_split_seq_params == "no"))){
    // Process Parameters for params.Parse_header_parse_headers:
    params.split_sequences_split_seq.field = "CONSCOUNT"
    params.split_sequences_split_seq.num = 20
}

if((params.edit_parse_header_table && (params.edit_parse_header_table == "no"))){
    // Process Parameters for params.Parse_header_table_parse_headers:
    params.Parse_header_table_parse_headers.method = "table"
    params.Parse_header_table_parse_headers.act = "min"
    params.Parse_header_table_parse_headers.args = "-f ID CREGION CONSCOUNT DUPCOUNT"
}


if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.R2_primers){params.R2_primers = ""} 
if (!params.R1_primers){params.R1_primers = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)

if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_2_reads_g0_0}
 } else {  
	g_2_reads_g0_0 = Channel.empty()
 }

Channel.value(params.mate).into{g_3_mate_g0_0;g_3_mate_g0_5;g_3_mate_g11_9;g_3_mate_g13_10;g_3_mate_g13_12;g_3_mate_g16_9;g_3_mate_g18_15;g_3_mate_g18_12;g_3_mate_g31_9;g_3_mate_g31_11}
g_6_R2_primers_g31_11 = params.R2_primers && file(params.R2_primers, type: 'any').exists() ? file(params.R2_primers, type: 'any') : ch_empty_file_2
g_7_R1_primers_g31_11 = params.R1_primers && file(params.R1_primers, type: 'any').exists() ? file(params.R1_primers, type: 'any') : ch_empty_file_1


process Filter_Sequence_Quality_filter_seq_quality {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /FS_.*$/) "FS_logs/$filename"}
input:
 set val(name),file(reads) from g_2_reads_g0_0
 val mate from g_3_mate_g0_0

output:
 set val(name), file("*_${method}-pass.fastq")  into g0_0_reads0_g31_11
 file "FS_*"  into g0_0_logFile1_g0_5

script:
method = params.Filter_Sequence_Quality_filter_seq_quality.method
nproc = params.Filter_Sequence_Quality_filter_seq_quality.nproc
q = params.Filter_Sequence_Quality_filter_seq_quality.q
n_length = params.Filter_Sequence_Quality_filter_seq_quality.n_length
n_missing = params.Filter_Sequence_Quality_filter_seq_quality.n_missing
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="quality"){
	q = "-q ${q}"
	n_length = ""
	n_missing = ""
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
	}else{
		q = ""
		n_length = ""
		n_missing = "-n ${n_missing}"
	}
}

readArray = reads.toString().split(' ')	

if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --outname ${name}_R1 --nproc ${nproc} --log FS_R1_${name}.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --outname ${name}_R2 --nproc ${nproc} --log FS_R2_${name}.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --outname ${name} --nproc ${nproc} --log FS_${name}.log
	"""
}


}


process Filter_Sequence_Quality_parse_log_FS {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "FS_log_tables/$filename"}
input:
 file log_file from g0_0_logFile1_g0_5
 val mate from g_3_mate_g0_5

output:
 file "*.tab"  into g0_5_logFile00

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f ID QUALITY
"""

}



process Mask_Primer_MaskPrimers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-pass.fastq$/) "MP_outputparam/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-fail.fastq$/) "MP_failed_reads/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /MP_.*$/) "MP_logs/$filename"}
input:
 val mate from g_3_mate_g31_11
 set val(name),file(reads) from g0_0_reads0_g31_11
 file R1_primers from g_7_R1_primers_g31_11
 file R2_primers from g_6_R2_primers_g31_11

output:
 set val(name), file("*_primers-pass.fastq")  into g31_11_reads0_g11_9
 set val(name), file("*_primers-fail.fastq") optional true  into g31_11_reads_failed11
 file "MP_*"  into g31_11_logFile2_g31_9

script:
method = params.Mask_Primer_MaskPrimers.method
barcode_field = params.Mask_Primer_MaskPrimers.barcode_field
primer_field = params.Mask_Primer_MaskPrimers.primer_field
barcode = params.Mask_Primer_MaskPrimers.barcode
revpr = params.Mask_Primer_MaskPrimers.revpr
mode = params.Mask_Primer_MaskPrimers.mode
failed = params.Mask_Primer_MaskPrimers.failed
nproc = params.Mask_Primer_MaskPrimers.nproc
maxerror = params.Mask_Primer_MaskPrimers.maxerror
umi_length = params.Mask_Primer_MaskPrimers.umi_length
start = params.Mask_Primer_MaskPrimers.start
extract_length = params.Mask_Primer_MaskPrimers.extract_length
maxlen = params.Mask_Primer_MaskPrimers.maxlen
skiprc = params.Mask_Primer_MaskPrimers.skiprc
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""

def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bf} ${pf} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk}")
}}

if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
    readArray = reads.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	"""
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed}
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed}
	"""
}else{
	args_1 = args_values[0]
	println args_1
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed}
	"""
}

}


process Pair_Sequence_pre_consensus_pair_seq {

input:
 set val(name),file(reads) from g31_11_reads0_g11_9
 val mate from g_3_mate_g11_9

output:
 set val(name),file("*_pair-pass.fastq")  into g11_9_reads0_g13_10

script:
coord = params.Pair_Sequence_pre_consensus_pair_seq.coord
act = params.Pair_Sequence_pre_consensus_pair_seq.act
copy_fields_1 = params.Pair_Sequence_pre_consensus_pair_seq.copy_fields_1
copy_fields_2 = params.Pair_Sequence_pre_consensus_pair_seq.copy_fields_2
failed = params.Pair_Sequence_pre_consensus_pair_seq.failed
nproc = params.Pair_Sequence_pre_consensus_pair_seq.nproc

if(mate=="pair"){
	
	act = (act=="none") ? "" : "--act ${act}"
	failed = (failed=="true") ? "--failed" : "" 
	copy_fields_1 = (copy_fields_1=="") ? "" : "--1f ${copy_fields_1}" 
	copy_fields_2 = (copy_fields_2=="") ? "" : "--2f ${copy_fields_2}" 
	
	readArray = reads.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	"""
	PairSeq.py -1 ${R1} -2 ${R2} ${copy_fields_1} ${copy_fields_2} --coord ${coord} ${act} ${failed}
	"""
}else{
	
	"""
	echo -e 'PairSeq works only on pair-end reads.'
	"""
}


}

boolean isCollectionOrArray_bc(object) {    
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}

def args_creator_bc(barcode_field, primer_field, act, copy_field, mincount, minqual, minfreq, maxerror, prcons, maxgap, maxdiv, dep){
	def args_values;
    if(isCollectionOrArray_bc(barcode_field) || isCollectionOrArray_bc(primer_field) || isCollectionOrArray_bc(copy_field) || isCollectionOrArray_bc(mincount) || isCollectionOrArray_bc(minqual) || isCollectionOrArray_bc(minfreq) || isCollectionOrArray_bc(maxerror) || isCollectionOrArray_bc(prcons) || isCollectionOrArray_bc(maxgap) || isCollectionOrArray_bc(maxdiv) || isCollectionOrArray_bc(dep)){
    	primer_field = (isCollectionOrArray_bc(primer_field)) ? primer_field : [primer_field,primer_field]
    	act = (isCollectionOrArray_bc(act)) ? act : [act,act]
    	copy_field = (isCollectionOrArray_bc(copy_field)) ? copy_field : [copy_field,copy_field]
    	mincount = (isCollectionOrArray_bc(mincount)) ? mincount : [mincount,mincount]
    	minqual = (isCollectionOrArray_bc(minqual)) ? minqual : [minqual,minqual]
    	minfreq = (isCollectionOrArray_bc(minfreq)) ? minfreq : [minfreq,minfreq]
    	maxerror = (isCollectionOrArray_bc(maxerror)) ? maxerror : [maxerror,maxerror]
    	prcons = (isCollectionOrArray_bc(prcons)) ? prcons : [prcons,prcons]
    	maxgap = (isCollectionOrArray_bc(maxgap)) ? maxgap : [maxgap,maxgap]
    	maxdiv = (isCollectionOrArray_bc(maxdiv)) ? maxdiv : [maxdiv,maxdiv]
    	dep = (isCollectionOrArray_bc(dep)) ? dep : [dep,dep]
    	args_values = []
        [barcode_field,primer_field,act,copy_field,mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep].transpose().each { bf,pf,a,cf,mc,mq,mf,mr,pc,mg,md,d -> {
            bf = (bf=="") ? "" : "--bf ${bf}"
            pf = (pf=="") ? "" : "--pf ${pf}" 
            a = (a=="none") ? "" : "--act ${a}" 
            cf = (cf=="") ? "" : "--cf ${cf}" 
            mr = (mr=="none") ? "" : "--maxerror ${mr}" 
            pc = (pc=="none") ? "" : "--prcons ${pc}" 
            mg = (mg=="none") ? "" : "--maxgap ${mg}" 
            md = (md=="none") ? "" : "--maxdiv ${md}" 
            d = (d=="true") ? "--dep" : "" 
            args_values.add("${bf} ${pf} ${a} ${cf} -n ${mc} -q ${mq} --freq ${mf} ${mr} ${pc} ${mg} ${md} ${d}")
        }}
    }else{
        barcode_field = (barcode_field=="") ? "" : "--bf ${barcode_field}"
        primer_field = (primer_field=="") ? "" : "--pf ${primer_field}" 
        act = (act=="none") ? "" : "--act ${act}" 
        copy_field = (copy_field=="") ? "" : "--cf ${copy_field}" 
        maxerror = (maxerror=="none") ? "" : "--maxerror ${maxerror}" 
        prcons = (prcons=="none") ? "" : "--prcons ${prcons}" 
        maxgap = (maxgap=="none") ? "" : "--maxgap ${maxgap}" 
        maxdiv = (maxdiv=="none") ? "" : "--maxdiv ${maxdiv}" 
        dep = (dep=="true") ? "--dep" : "" 
        args_values = "${barcode_field} ${primer_field} ${act} ${copy_field} -n ${mincount} -q ${minqual} --freq ${minfreq} ${maxerror} ${prcons} ${maxgap} ${maxdiv} ${dep}"
    }
    return args_values
}


process Build_Consensus_build_consensus {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_consensus-pass.fastq$/) "BC_pass_read/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /BC_.*$/) "BC_log/$filename"}
input:
 set val(name),file(reads) from g11_9_reads0_g13_10
 val mate from g_3_mate_g13_10

output:
 set val(name),file("*_consensus-pass.fastq")  into g13_10_reads0_g16_9
 file "BC_*"  into g13_10_logFile1_g13_12

script:
failed = params.Build_Consensus_build_consensus.failed
nproc = params.Build_Consensus_build_consensus.nproc
barcode_field = params.Build_Consensus_build_consensus.barcode_field
primer_field = params.Build_Consensus_build_consensus.primer_field
act = params.Build_Consensus_build_consensus.act
copy_field = params.Build_Consensus_build_consensus.copy_field
mincount = params.Build_Consensus_build_consensus.mincount
minqual = params.Build_Consensus_build_consensus.minqual
minfreq = params.Build_Consensus_build_consensus.minfreq
maxerror = params.Build_Consensus_build_consensus.maxerror
prcons = params.Build_Consensus_build_consensus.prcons
maxgap = params.Build_Consensus_build_consensus.maxgap
maxdiv = params.Build_Consensus_build_consensus.maxdiv
dep = params.Build_Consensus_build_consensus.dep
//* @style @condition:{act="none",},{act="min",copy_field},{act="max",copy_field},{act="sum",copy_field},{act="set",copy_field},{act="majority",copy_field} @array:{barcode_field,primer_field,act,copy_field,mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep} @multicolumn:{failed,nproc},{barcode_field,primer_field,act,copy_field}, {mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep}

args_values_bc = args_creator_bc(barcode_field, primer_field, act, copy_field, mincount, minqual, minfreq, maxerror, prcons, maxgap, maxdiv, dep)

// args 
if(isCollectionOrArray_bc(args_values_bc)){
	args_1 = args_values_bc[0]
	args_2 = args_values_bc[1]
}else{
	args_1 = args_values_bc
	args_2 = args_values_bc
}


failed = (failed=="true") ? "--failed" : "" 
if(mate=="pair"){
	// files
	readArray = reads.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	"""
	BuildConsensus.py -s $R1 ${args_1} --outname ${name}_R1 --log BC_${name}_R1.log ${failed} --nproc ${nproc}
	BuildConsensus.py -s $R2 ${args_2} --outname ${name}_R2 --log BC_${name}_R2.log ${failed} --nproc ${nproc}
	"""
}else{
	"""
	BuildConsensus.py -s $reads ${args_1} --outname ${name} --log BC_${name}.log ${failed} --nproc ${nproc}
	"""
}


}


process Build_Consensus_parse_log_BC_copy {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "BC_log_table/$filename"}
input:
 file log_file from g13_10_logFile1_g13_12
 val mate from g_3_mate_g13_12

output:
 file "*.tab"  into g13_12_logFile00

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f BARCODE SEQCOUNT CONSCOUNT PRCONS PRFREQ ERROR
"""

}


process Pair_Sequence_post_consensus_pair_seq {

input:
 set val(name),file(reads) from g13_10_reads0_g16_9
 val mate from g_3_mate_g16_9

output:
 set val(name),file("*_pair-pass.fastq")  into g16_9_reads0_g18_12

script:
coord = params.Pair_Sequence_post_consensus_pair_seq.coord
act = params.Pair_Sequence_post_consensus_pair_seq.act
copy_fields_1 = params.Pair_Sequence_post_consensus_pair_seq.copy_fields_1
copy_fields_2 = params.Pair_Sequence_post_consensus_pair_seq.copy_fields_2
failed = params.Pair_Sequence_post_consensus_pair_seq.failed
nproc = params.Pair_Sequence_post_consensus_pair_seq.nproc

if(mate=="pair"){
	
	act = (act=="none") ? "" : "--act ${act}"
	failed = (failed=="true") ? "--failed" : "" 
	copy_fields_1 = (copy_fields_1=="") ? "" : "--1f ${copy_fields_1}" 
	copy_fields_2 = (copy_fields_2=="") ? "" : "--2f ${copy_fields_2}" 
	
	readArray = reads.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	"""
	PairSeq.py -1 ${R1} -2 ${R2} ${copy_fields_1} ${copy_fields_2} --coord ${coord} ${act} ${failed}
	"""
}else{
	
	"""
	echo -e 'PairSeq works only on pair-end reads.'
	"""
}


}


process Assemble_pairs_assemble_pairs {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_assemble-pass.f.*$/) "AP_pass_read/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /AP_.*$/) "AP_log_files/$filename"}
input:
 set val(name),file(reads) from g16_9_reads0_g18_12
 val mate from g_3_mate_g18_12

output:
 set val(name),file("*_assemble-pass.f*")  into g18_12_reads0_g21_15
 file "AP_*"  into g18_12_logFile1_g18_15
 set val(name),file("*_assemble-fail.f*") optional true  into g18_12_reads_failed22

script:
method = params.Assemble_pairs_assemble_pairs.method
coord = params.Assemble_pairs_assemble_pairs.coord
rc = params.Assemble_pairs_assemble_pairs.rc
head_fields_R1 = params.Assemble_pairs_assemble_pairs.head_fields_R1
head_fields_R2 = params.Assemble_pairs_assemble_pairs.head_fields_R2
failed = params.Assemble_pairs_assemble_pairs.failed
fasta = params.Assemble_pairs_assemble_pairs.fasta
nproc = params.Assemble_pairs_assemble_pairs.nproc
alpha = params.Assemble_pairs_assemble_pairs.alpha
maxerror = params.Assemble_pairs_assemble_pairs.maxerror
minlen = params.Assemble_pairs_assemble_pairs.minlen
maxlen = params.Assemble_pairs_assemble_pairs.maxlen
scanrev = params.Assemble_pairs_assemble_pairs.scanrev
minident = params.Assemble_pairs_assemble_pairs.minident
evalue = params.Assemble_pairs_assemble_pairs.evalue
maxhits = params.Assemble_pairs_assemble_pairs.maxhits
fill = params.Assemble_pairs_assemble_pairs.fill
aligner = params.Assemble_pairs_assemble_pairs.aligner
// align_exec = params.Assemble_pairs_assemble_pairs.// align_exec
// dbexec = params.Assemble_pairs_assemble_pairs.// dbexec
gap = params.Assemble_pairs_assemble_pairs.gap
usearch_version = params.Assemble_pairs_assemble_pairs.usearch_version
//* @style @condition:{method="align",alpha,maxerror,minlen,maxlen,scanrev}, {method="sequential",alpha,maxerror,minlen,maxlen,scanrev,ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="reference",ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="join",gap} @multicolumn:{method,coord,rc,head_fields_R1,head_fields_R2,failed,nrpoc,usearch_version},{alpha,maxerror,minlen,maxlen,scanrev}, {ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec}, {gap} 

// args
coord = "--coord ${coord}"
rc = "--rc ${rc}"
head_fields_R1 = (head_fields_R1!="") ? "--1f ${head_fields_R1}" : ""
head_fields_R2 = (head_fields_R2!="") ? "--2f ${head_fields_R2}" : ""
failed = (failed=="false") ? "" : "--failed"
fasta = (fasta=="false") ? "" : "--fasta"
nproc = "--nproc ${nproc}"

scanrev = (scanrev=="false") ? "" : "--scanrev"
fill = (fill=="false") ? "" : "--fill"

// align_exec = (align_exec!="") ? "--exec ${align_exec}" : ""
// dbexec = (dbexec!="") ? "--dbexec ${dbexec}" : ""

def ref_file;
if(params.assemble_ref != null) ref_file = (assemble_reference.startsWith('NO_FILE')) ? "" : "-r ${assemble_reference}"

args = ""

if(method=="align"){
	args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev}"
}else{
	if(method=="sequential"){
		args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev} ${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
	}else{
		if(method=="reference"){
			args = "${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
		}else{
			args = "--gap ${gap}"
		}
	}
}


readArray = reads.toString().split(' ')	
if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	"""
	if [ "${method}" != "align" ]; then
		if  [ "${aligner}" == "usearch" ]; then
			wget -q --show-progress --no-check-certificate https://drive5.com/downloads/usearch${usearch_version}_i86linux32.gz
			gunzip usearch${usearch_version}_i86linux32.gz
			chmod +x usearch${usearch_version}_i86linux32
			mv usearch${usearch_version}_i86linux32 /usr/local/bin/usearch2
			align_exec="--exec /usr/local/bin/usearch2"
			dbexec="--dbexec /usr/local/bin/usearch2"
		else
			align_exec="--exec /usr/local/bin/blastn"
			dbexec="--dbexec /usr/local/bin/makeblastdb"
		fi
	else
		align_exec=""
		dbexec=""
	fi
	
	echo \$align_exec
	echo \$dbexec
	
	
	
	AssemblePairs.py ${method} -1 ${R1} -2 ${R2} ${coord} ${rc} ${head_fields_R1} ${head_fields_R2} ${args} \$align_exec \$dbexec ${fasta} ${failed} --log AP_${name}.log ${nproc}
	"""

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}

}


process Assemble_pairs_parse_log_AP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "AP_log_table/$filename"}
input:
 file log_file from g18_12_logFile1_g18_15
 val mate from g_3_mate_g18_15

output:
 file "*.tab"  into g18_15_logFile00

script:
field_to_parse = params.Assemble_pairs_parse_log_AP.field_to_parse
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ${field_to_parse}
"""


}


process Parse_header_parse_headers {

input:
 set val(name), file(reads) from g18_12_reads0_g21_15

output:
 set val(name),file("*${out}")  into g21_15_reads0_g22_16

script:
method = params.Parse_header_parse_headers.method
act = params.Parse_header_parse_headers.act
args = params.Parse_header_parse_headers.args

out="_reheader.fastq"
if(method=="collapse" || method=="add" || method=="copy" || method=="rename" || method=="merge"){
	"""
	ParseHeaders.py  ${method} -s ${reads} ${args} --act ${act}
	"""
}else{
	if(method=="table"){
			out=".tab"
			"""
			ParseHeaders.py ${method} -s ${reads} -o ${name}.tab ${args}
			"""	
	}else{
		"""
		ParseHeaders.py ${method} -s ${reads} ${args}
		"""		
	}
}



}


process collapse_sequences_collapse_seq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_collapse-unique.fast.*$/) "reads_unique/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_collapse-duplicate.fast.*$/) "reads_duplicated/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_collapse-undetermined.fast.*$/) "reads_unique/$filename"}
input:
 set val(name), file(reads) from g21_15_reads0_g22_16

output:
 set val(name),  file("*_collapse-unique.fast*")  into g22_16_reads0_g26_20
 set val(name),  file("*_collapse-duplicate.fast*") optional true  into g22_16_reads_duplicate11
 set val(name),  file("*_collapse-undetermined.fast*") optional true  into g22_16_reads_undetermined22

script:
max_missing = params.collapse_sequences_collapse_seq.max_missing
inner = params.collapse_sequences_collapse_seq.inner
fasta = params.collapse_sequences_collapse_seq.fasta
act = params.collapse_sequences_collapse_seq.act
uf = params.collapse_sequences_collapse_seq.uf
cf = params.collapse_sequences_collapse_seq.cf
nproc = params.collapse_sequences_collapse_seq.nproc

inner = (inner=="true") ? "--inner" : ""
fasta = (fasta=="true") ? "--fasta" : ""
act = (act=="none") ? "" : "--act ${act}"
cf = (cf=="") ? "" : "--cf ${cf}"
uf = (uf=="") ? "" : "--uf ${uf}"

"""
CollapseSeq.py -s ${reads} -n ${max_missing} ${fasta} ${inner} ${uf} ${cf} ${act}
"""

}


process split_sequences_split_seq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_atleast-.*.fastq$/) "split_sequence_reads/$filename"}
input:
 set val(name),file(reads) from g22_16_reads0_g26_20

output:
 set val(name), file("*_atleast-*.fastq")  into g26_20_reads0_g28_15

script:
field = params.split_sequences_split_seq.field
num = params.split_sequences_split_seq.num

if(num!=0){
	num = " --num ${num}"
}else{
	num = ""
}

"""
SplitSeq.py group -s ${reads} -f ${field} ${num}
"""

}


process Parse_header_table_parse_headers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*${out}$/) "parse_header_table/$filename"}
input:
 set val(name), file(reads) from g26_20_reads0_g28_15

output:
 set val(name),file("*${out}")  into g28_15_reads00

script:
method = params.Parse_header_table_parse_headers.method
act = params.Parse_header_table_parse_headers.act
args = params.Parse_header_table_parse_headers.args

out="_reheader.fastq"
if(method=="collapse" || method=="add" || method=="copy" || method=="rename" || method=="merge"){
	"""
	ParseHeaders.py  ${method} -s ${reads} ${args} --act ${act}
	"""
}else{
	if(method=="table"){
			out=".tab"
			"""
			ParseHeaders.py ${method} -s ${reads} -o ${name}.tab ${args}
			"""	
	}else{
		"""
		ParseHeaders.py ${method} -s ${reads} ${args}
		"""		
	}
}



}


process Mask_Primer_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "MP_logs/$filename"}
input:
 val mate from g_3_mate_g31_9
 file log_file from g31_11_logFile2_g31_9

output:
 file "*.tab"  into g31_9_logFile00

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
