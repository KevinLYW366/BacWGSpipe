####################################################################
# Script to extract and merge results for all samples              #
# This script is designed for different analysis steps in workflow #
# By Yewei on 2022.11.4                                            #
####################################################################

def get_options():
    import argparse
    # create the top-level parser
    description = "Extract and merge results for all samples"
    parser = argparse.ArgumentParser(description=description,
                                     prog='merge_results.py')
    parser.add_argument('-m', '--module', default=False,
                        choices=['rgi', 'digis', 'icefinder', 'islandpath',
                                 'phigaro', 'plasmidfinder', 'platon', 'vfdb',
                                 'mob_suite', 'plsdb', 'mlst', 'cctyper',
                                 'abricate_card'],
                        help='module with results to extract and merge')
    parser.add_argument('-i', '--result-dir', default=False,
                        help='path to result directory, please match this path with the analysis module input')
    parser.add_argument('-o', '--output-dir', default=False,
                        help='path to output directory under which merged results file will be generated')
    parser.add_argument('-s', '--samplelist', default=False,
                        help='sample list file')
    return parser.parse_args()


# Function to extract pos in genome (contig, start, end, orientation) from gff file based on gene_id
def convert_gene_id(gene_id, gff):
    pattern = "ID={};".format(gene_id)
    gff_lines = [line.strip() for line in open(gff)]
    id_line = list(filter(lambda x: pattern in x, gff_lines))[0]
    contig = id_line.split('\t')[0]
    start = id_line.split('\t')[3]
    end = id_line.split('\t')[4]
    orientation = id_line.split('\t')[6]
    return contig, start, end, orientation


# Function to read sample list
def read_sample_list(sample_list_file):
    sample_list = [line.strip() for line in open(sample_list_file)]
    return sample_list


# Function to extract digIS result and output
def digis_results(output_dir, result_dir, sample_list):
    result_output = os.path.join(output_dir, "all_digis_results.xls")
    with open(result_output, 'w') as f_out:
        f_out.write(
            "sample\tpos\tid\tlevel\tqid\tqstart\tqend\tsid\tsstart\tsend\tstrand\tacc\tscore\tevalue\tORF_sim\tIS_sim\tGenBank_class\n")
        for sample in sample_list:
            digis_result = os.path.join(result_dir, "{}/digis/results/{}.csv".format(sample, sample))
            with open(digis_result, 'r') as f_in:
                f_in.readline()
                for line in f_in:
                    content = line.split(',')
                    contig, start, end = content[5], content[6], content[7]
                    pos = "{}:{}-{}".format(contig, start, end)
                    f_out.write("{}\t{}\t{}".format(sample, pos, line.replace(',', '\t')))


# Function to extract ICEfinder result and output
def icefinder_results(output_dir, result_dir, sample_list):
    result_output = os.path.join(output_dir, "all_icefinder_results.xls")
    with open(result_output, 'w') as f_out:
        f_out.write(
            "sample\tpos\tICEfinder_Job_id\tStrain\tGenome_len\tICEfinder_output_name\tDescription\tCoordinate\tLength\toriT\tGC\tGenome_GC\t|Delta_GC|\tARG\tVF\n")
        for sample in sample_list:
            icefinder_result = os.path.join(result_dir, "{}/icefinder/{}_icefinder_result.txt".format(sample, sample))
            with open(icefinder_result, 'r') as f_in:
                f_in.readline()
                for line in f_in:
                    content = line.split('\t')
                    contig = content[0]
                    start, end = content[5].split('..')
                    pos = "{}:{}-{}".format(contig, start, end)
                    f_out.write("{}\t{}\t{}".format(sample, pos, line))


# Function to extract IslandPath-DIMOB result and output
def islandpath_results(output_dir, result_dir, sample_list):
    result_output = os.path.join(output_dir, "all_islandpath_results.xls")
    with open(result_output, 'w') as f_out:
        f_out.write("sample\tpos\tcontig\tstart\tend\n")
        for sample in sample_list:
            islandpath_result = os.path.join(result_dir, "{}/islandpath/{}_GIs.bed".format(sample, sample))
            with open(islandpath_result, 'r') as f_in:
                for line in f_in:
                    content = line.strip().split('\t')
                    contig, start, end = content[0], content[1], content[2]
                    pos = "{}:{}-{}".format(contig, start, end)
                    f_out.write("{}\t{}\t{}".format(sample, pos, line))


# Function to extract Phigaro result and output
def phigaro_results(output_dir, result_dir, sample_list):
    result_output = os.path.join(output_dir, "all_phigaro_results.xls")
    with open(result_output, 'w') as f_out:
        f_out.write("sample\tpos\tscaffold\tbegin\tend\ttransposable\ttaxonomy\tvog\n")
        for sample in sample_list:
            phigaro_result = os.path.join(result_dir, "{}/phigaro/{}.phigaro.tsv".format(sample, sample))
            with open(phigaro_result, 'r') as f_in:
                f_in.readline()
                for line in f_in:
                    content = line.split('\t')
                    pos = "{}:{}-{}".format(content[0], content[1], content[2])
                    f_out.write("{}\t{}\t{}".format(sample, pos, line))


# Function to extract Plasmidfinder result and output
def plasmidfinder_results(output_dir, result_dir, sample_list):
    result_output = os.path.join(output_dir, "all_plasmidfinder_results.xls")
    with open(result_output, 'w') as f_out:
        f_out.write(
            "sample\tpos\tDatabase\tPlasmid\tIdentity\tQuery / Template length\tContig\tPosition in contig\tNote\tAccession number\n")
        for sample in sample_list:
            plasmidfinder_result = os.path.join(result_dir, "{}/plasmidfinder/results_tab.tsv".format(sample))
            with open(plasmidfinder_result, 'r') as f_in:
                f_in.readline()
                for line in f_in:
                    content = line.split('\t')
                    contig = content[4]
                    start, end = content[5].split('..')
                    pos = "{}:{}-{}".format(contig, start, end)
                    f_out.write("{}\t{}\t{}".format(sample, pos, line))


# Function to extract RGI result and output
def rgi_results(output_dir, result_dir, sample_list):
    result_output = os.path.join(output_dir, "all_rgi_results.xls")
    with open(result_output, 'w') as f_out:
        f_out.write(
            "sample\tORF_ID\tReplicon\tStart\tStop\tOrientation\tCut_Off\tPass_Bitscore\tBest_Hit_Bitscore\tBest_Hit_ARO\tBest_Identities\tARO\tModel_type\tSNPs_in_Best_Hit_ARO\tOther_SNPs\tDrug Class\tResistance Mechanism\tAMR Gene Family\tPredicted_DNA\tPredicted_Protein	CARD_Protein_Sequence\tPercentage Length of Reference Sequence\tID\tModel_ID\tNudged\tNote\n")
        for sample in sample_list:
            gff = "results/03.genome_component/{}/prokka/{}.gff".format(sample, sample)
            rgi_result = os.path.join(result_dir, "{}/CARD/rgi/{}.report.txt".format(sample, sample))
            with open(rgi_result, 'r') as f_in:
                f_in.readline()
                for line in f_in:
                    items = line.split("\t")
                    gene_id = items[0].split(" ")[0]
                    contig, start, end, orientation = convert_gene_id(gene_id, gff)
                    items[1] = contig
                    items[2] = start
                    items[3] = end
                    items[4] = orientation
                    f_out.write("{}\t{}".format(sample, "\t".join(items)))


# Function to extract Platon result and output
def platon_results(output_dir, result_dir, sample_list):
    result_output = os.path.join(output_dir, "all_platon_results.xls")
    with open(result_output, 'w') as f_out:
        f_out.write(
            "sample\tpos\tID\tLength\tCoverage\t# ORFs\tRDS\tCircular\tInc Type(s)\t# Replication\t# Mobilization\t# OriT\t# Conjugation\t# AMRs\t# rRNAs\t# Plasmid Hits\n")
        for sample in sample_list:
            platon_result = os.path.join(result_dir, "{}/platon/{}.tsv".format(sample, sample))
            with open(platon_result, 'r') as f_in:
                f_in.readline()
                for line in f_in:
                    content = line.split('\t')
                    contig = content[0]
                    # Platon identify the complete contig as plasmid
                    start = 1
                    end = content[1]
                    pos = "{}:{}-{}".format(contig, start, end)
                    f_out.write("{}\t{}\t{}".format(sample, pos, line))


# Function to extract VFDB result and output
def vfdb_results(output_dir, result_dir, sample_list):
    result_output = os.path.join(output_dir, "all_vfdb_results.xls")
    with open(result_output, 'w') as f_out:
        f_out.write(
            "sample\tqseqid\treplicon\tstart\tend\torientation\tqstart\tqend\tqlen\tsseqid\tsstart\tsend\tslen\tevalue\tlength\tpident\tgaps\tgapopen\tstitle\n")
        for sample in sample_list:
            gff = "results/03.genome_component/{}/prokka/{}.gff".format(sample, sample)
            vfdb_result = os.path.join(result_dir, "{}/VFDB/{}_vfdb_blastn_onGenes.txt".format(sample, sample))
            with open(vfdb_result, 'r') as f_in:
                f_in.readline()
                for line in f_in:
                    items = line.split("\t")
                    gene_id = items[0]
                    contig, start, end, orientation = convert_gene_id(gene_id, gff)
                    line = "\t".join([gene_id, contig, start, end, orientation] + items[1:])
                    f_out.write("{}\t{}".format(sample, line))


# Function to extract MOB-suite result and output
def mob_suite_results(output_dir, result_dir, sample_list):
    result_output = os.path.join(output_dir, "all_mob-suite_results.xls")
    with open(result_output, 'w') as f_out:
        n = 1
        for sample in sample_list:
            mob_suite_result = os.path.join(result_dir, "{}/mob-suite/{}_mobtyper_result.xls".format(sample, sample))
            with open(mob_suite_result, 'r') as f_in:
                header = f_in.readline()
                if n:
                    header_new = ["sample_id", "plasmid_id"] + header.split("\t")[1:]
                    f_out.write("\t".join(header_new))
                    n = 0
                for line in f_in:
                    f_out.write("{}\t{}".format(sample, line))


# Function to extract PLSDB annotation result and output
def plsdb_results(output_dir, result_dir, sample_list):
    result_output = os.path.join(output_dir, "all_plsdb_results.xls")
    with open(result_output, 'w') as f_out:
        header = ""
        for sample in sample_list:
            sample_dir = os.path.join(result_dir, "{}/plsdb".format(sample))
            filename_list = os.listdir(sample_dir)
            filename_list.sort()
            for filename in filename_list:
                if "plas" in filename:
                    with open(os.path.join(sample_dir, filename), 'r') as f_in:
                        if header == "":
                            header = f_in.readline()
                            f_out.write(header)
                        else:
                            f_in.readline()
                        for line in f_in:
                            f_out.write(line)


# Function to extract MLST result and output
def mlst_results(output_dir, result_dir, sample_list):
    result_output = os.path.join(output_dir, "all_mlst_results.xls")
    with open(result_output, 'w') as f_out:
        header = "Sample\tPubMLST_scheme\tST\tAllele_IDs\n"
        f_out.write(header)
        for sample in sample_list:
            mlst_result = os.path.join(result_dir, "{}/mlst/{}_mlst.txt".format(sample, sample))
            with open(mlst_result, 'r') as f_in:
                line = f_in.readline()
                items = line.strip().split("\t")
                scheme = items[1]
                st = items[2]
                alleles = "; ".join(items[3:])
                f_out.write("{}\t{}\t{}\t{}\n".format(sample, scheme, st, alleles))


# Function to extract cctyper result and output
def cctyper_results(output_dir, result_dir, sample_list):
    result_output = os.path.join(output_dir, "all_cctyper_results.xls")
    with open(result_output, 'w') as f_out:
        f_out.write("Sample\tContig\tCRISPR\tStart\tEnd\tConsensus_repeat\tN_repeats\tRepeat_len\tSpacer_len_avg\tRepeat_identity\tSpacer_identity\tSpacer_len_sem\tTrusted\tPrediction\tSubtype\tSubtype_probability\n")
        for sample in sample_list:
            cctyper_result = os.path.join(result_dir, "{}/cctyper/results/crisprs_all.tab".format(sample))
            if os.path.isfile(cctyper_result):
                with open(cctyper_result, 'r') as f_in:
                    f_in.readline()
                    f_in.readline()
                    for line in f_in:
                        items = line.split("\t")
                        line = "\t".join([sample] + items)
                        f_out.write(line)


# Function to extract Abricate (CARD) result and output
def abricate_card_results(output_dir, result_dir, sample_list):
    result_output = os.path.join(output_dir, "all_abricate_card_results.xls")
    with open(result_output, 'w') as f_out:
        f_out.write(
            "SAMPLE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n")
        for sample in sample_list:
            abricate_card_result = os.path.join(result_dir, "{}/CARD/abricate/{}_abricate_card_results.xls".format(sample, sample))
            with open(abricate_card_result, 'r') as f_in:
                f_in.readline()
                for line in f_in:
                    items = line.split("\t")
                    f_out.write("{}\t{}".format(sample, "\t".join(items[1:])))


if __name__ == "__main__":
    # import library
    import sys
    import os

    # parse options
    options = get_options()

    # Read sample list
    sample_list = read_sample_list(options.samplelist)

    # extract and merge results for the input module
    # Available modules include rgi, digis, icefinder, islandpath,
    #   phigaro, plasmidfinder, platon, vfdb, mob_suite, plsdb
    if options.module == "rgi":
        rgi_results(options.output_dir, options.result_dir, sample_list)
    elif options.module == "digis":
        digis_results(options.output_dir, options.result_dir, sample_list)
    elif options.module == "icefinder":
        icefinder_results(options.output_dir, options.result_dir, sample_list)
    elif options.module == "islandpath":
        islandpath_results(options.output_dir, options.result_dir, sample_list)
    elif options.module == "phigaro":
        phigaro_results(options.output_dir, options.result_dir, sample_list)
    elif options.module == "plasmidfinder":
        plasmidfinder_results(options.output_dir, options.result_dir, sample_list)
    elif options.module == "platon":
        platon_results(options.output_dir, options.result_dir, sample_list)
    elif options.module == "vfdb":
        vfdb_results(options.output_dir, options.result_dir, sample_list)
    elif options.module == "mob_suite":
        mob_suite_results(options.output_dir, options.result_dir, sample_list)
    elif options.module == "plsdb":
        plsdb_results(options.output_dir, options.result_dir, sample_list)
    elif options.module == "mlst":
        mlst_results(options.output_dir, options.result_dir, sample_list)
    elif options.module == "cctyper":
        cctyper_results(options.output_dir, options.result_dir, sample_list)
    elif options.module == "abricate_card":
        abricate_card_results(options.output_dir, options.result_dir, sample_list)
