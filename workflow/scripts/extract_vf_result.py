############################################
# Script to extract and merge VFDB results #
# By Yewei on 2022.3.10                    #
############################################

# Output format:
#          VF1               VF2  ...
# Sample1  Contig:Start:End  ...
# Sample2  ...
# ...

# Library
import sys
import os
import re
import pandas as pd
import xlsxwriter

# Main
if __name__ == "__main__":
    # Input
    # Merged RGI and Abricate (CARD) results file
    vfdb_results_file = sys.argv[1]

    # Output
    matrix_output_file = sys.argv[2]
    vf_info_output_file = sys.argv[3]
    # Highlighted matrix with 0/1
    binary_matrix_output_file = sys.argv[4]

    # Extract VFDB result
    result_dict = {}
    vf_info_dict = {}

    with open(vfdb_results_file, 'r') as fin:
        # header line
        buf = fin.readline()
        buf = fin.readline()
        while buf:
            items = buf.strip().split('\t')
            sample, gene_id, contig, start, end, vf_info = items[0], items[1], items[2], items[3], items[4], items[18]
            pos = "{}:{}-{}".format(contig, start, end)
            m = re.search(r'(.*?\(.*?\)) \((.*?)\) (.*?)$', vf_info)
            if not m:
                m = re.search(r'(.*?) \((.*?)\) (.*?)$', vf_info)
            vf_vfdb_id = m.group(1)
            vf_gene = m.group(2)
            vf_info = m.group(3)
            if sample not in result_dict.keys():
                result_dict[sample] = {}
            if vf_gene not in result_dict[sample].keys():
                result_dict[sample][vf_gene] = pos
            else:
                result_dict[sample][vf_gene] = "{},{}".format(result_dict[sample][vf_gene], contig)
            # extract ARO info
            if vf_vfdb_id not in vf_info_dict.keys():
                vf_info_dict[vf_vfdb_id] = [vf_gene, vf_info]
            buf = fin.readline()


    # Convert dict to dataframe and output
    result_output_df = pd.DataFrame.from_dict(result_dict, orient='index')
    ## the order of samples in output follows that in sample list
    #result_output_df = result_output_df.reindex(sample_list, axis=0)
    #result_output_df = result_output_df.sort_index(axis=1)
    result_output_df.to_csv(matrix_output_file, header=True, index=True, sep='\t')
    # Highlighted matrix with 0/1
    result_output_df_binary = result_output_df.notna().astype(int)
    # result_output_df_binary.to_csv(binary_matrix_output_file, header=True, index=True, sep='\t')
    writer = pd.ExcelWriter(binary_matrix_output_file, engine='xlsxwriter')
    result_output_df_binary.to_excel(writer, sheet_name="Sheet1")
    workbook = writer.book
    worksheet = writer.sheets["Sheet1"]
    len_col = len(result_output_df_binary.columns)
    len_row = len(result_output_df_binary.index)
    # Light red fill with dark red text.
    format1 = workbook.add_format({'bg_color': '#FFC7CE',
                                   'font_color': '#9C0006'})
    worksheet.conditional_format(0, 0, len_row, len_col, {'type': 'cell',
                                                          'criteria': '=',
                                                          'value': 1,
                                                          'format': format1})
    writer.save()

    vf_info_output_df = pd.DataFrame.from_dict(vf_info_dict,
                                                orient='index',
                                                columns=["VF_gene",
                                                         "VF_description"])
    vf_info_output_df = vf_info_output_df.sort_index(axis=0)
    vf_info_output_df.index.name = "VFDB_internal_ID"
    vf_info_output_df.to_csv(vf_info_output_file, header=True, index=True, sep='\t')
