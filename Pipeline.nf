process ref_vcf2chrom {
	//debug true

	executor "local"

	input:
	tuple file(vcf), file(vcf_index)

	output:
	tuple stdout, file(vcf), file(vcf_index)

	"""
	n_chrom=`bcftools index -s ${vcf} | wc -l`
	if [[ \${n_chrom} -gt 1 ]]; then
		echo "Multiple chromosomes within one reference panel VCF are not allowed." 1>&2
		exit 1
	fi
	chrom=`bcftools index -s ${vcf} | cut -f1`
	printf "\${chrom}"
	"""
}


process study_vcf2chrom {
	//debug true

	executor "local"

	input:
	tuple file(vcf), file(vcf_index)

	output:
	tuple stdout, file(vcf), file(vcf_index)

	"""
	n_chrom=`bcftools index -s ${vcf} | wc -l`
	if [[ \${n_chrom} -gt 1 ]]; then
		echo "Multiple chromosomes within one study VCF are not allowed." 1>&2
		exit 1
	fi
	chrom=`bcftools index -s ${vcf} | cut -f1`
	printf "\${chrom}"
	"""
}


process infer_auto_chrom {
	label "INFERENCE"

	//debug true
	//executor "local"

	cache "lenient"
        //scratch true
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return "retry" }
	maxRetries 0

	input:
	tuple val(chrom), file(study_vcf), file(study_vcf_index), file(ref_vcf), file(ref_vcf_index), file(genetic_map)
	each file(sample_map)

        output:
	file "${chrom}.rfmix_inference.rfmix.Q"
        file "${chrom}.rfmix_inference.fb.tsv"
	file "${chrom}.rfmix_inference.msp.tsv"
	file "${chrom}.rfmix_inference.sis.tsv"

	publishDir "Results", pattern: "${chrom}.rfmix_inference.*", mode: "copy"

	"""
	${params.rfmix2_cmd} -f ${study_vcf} -r ${ref_vcf} -g ${genetic_map} -m ${sample_map} --chromosome=${chrom} -o ${chrom}.rfmix_inference
	"""
}


process merge_q {
        //debug true

        executor "local"

	input:
        file q_files

	output:
	file("ALL.rfmix_inference.rfmix.Q")

	publishDir "Results", pattern: "ALL.rfmix_inference.rfmix.Q", mode: "copy"

	"""
	#!/usr/bin/env python3
	import pandas as pd
	import glob
	df_merged = pd.concat([ pd.read_csv(f, sep = '\t', skiprows = 1) for f in glob.glob('*.Q') ]).groupby(['#sample']).mean().reset_index()
	df_merged.to_csv('ALL.rfmix_inference.rfmix.Q', sep = '\t', index = False)
	"""
}


workflow {
	ref_vcfs = Channel.fromPath(params.reference_vcfs).map{ vcf -> [ vcf, vcf + (file(vcf + ".tbi").exists() ? ".tbi" : ".csi") ] }
        study_vcfs = Channel.fromPath(params.study_vcfs).map{ vcf -> [ vcf, vcf + (file(vcf +  ".tbi").exists() ? ".tbi" : ".csi") ] }
	sample_map = Channel.fromPath(params.reference_sample_map)

	ref_vcf2chrom(ref_vcfs)
	study_vcf2chrom(study_vcfs)

	// split autosomals and X chromosome
	ref_vcf2chrom.out.branch {
		auto_chroms: it[0] =~ /^(chr?)[1-9][0-9]*$/
		x_chrom: it[0] =~ /^(chr?)X$/
	}.set{ ref }

	study_vcf2chrom.out.branch {
		auto_chroms: it[0] =~ /^(chr?)[1-9][0-9]*$/
		x_chrom: it[0] =~ /^(chr?)X$/
	}.set{ study }

	// group  study and reference files by chromosome
	study_ref_auto_chroms = study.auto_chroms.combine(ref.auto_chroms, by: 0).map { it -> it + [ file("$workflow.projectDir/Genetic_maps/${it[0]}.GRCh38.map") ] }

	infer_auto_chrom(study_ref_auto_chroms, sample_map)
	
	merge_q(infer_auto_chrom.out[0].collect())
}

