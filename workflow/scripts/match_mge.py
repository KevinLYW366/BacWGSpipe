##################################################################
# Script to match AMR or VF with MGE in the same assembly contig #
# By Yewei on 2022.3.11                                          #
##################################################################

# Function to extract contig, start, end from pos (contig:start-end)
def get_contig(pos):
    return pos.split(':')[0]


def get_start(pos):
    start = int(pos.split(':')[1].split("-")[0])
    print("{};{}".format(pos, start))
    return start


def get_end(pos):
    end = int(pos.split(':')[1].split("-")[1])
    print("{};{}".format(pos, end))
    return end


# Function to match MGE results
def match_mge(sample, pos, result_file, distance):
    mge_pos_list = []
    with open(result_file, 'r') as f:
        f.readline()
        for p in pos.split(","):
            for line in f:
                content = line.split('\t')
                if sample == content[0] and get_contig(p) == get_contig(content[1]):
                    if distance >= 0:
                        if not ((get_start(p) - get_end(content[1]) > distance)
                                or (get_start(content[1]) - get_end(p) > distance)):
                            mge_pos_list.append(content[1])
                    elif distance == -1:
                        mge_pos_list.append(content[1])
    return mge_pos_list


# Main
if __name__ == "__main__":
    # Library
    import sys
    import os
    import re
    import pandas as pd
    import math

    # Input
    # AMR matrix for all samples
    amr_matrix_file = sys.argv[1]
    # VF matrix for all samples
    vf_matrix_file = sys.argv[2]
    # Association between AMR and MGE results output
    amr_mge_output = sys.argv[3]
    # Association between VF and MGE results output
    vf_mge_output = sys.argv[4]
    # Merged phigaro results
    phigaro_results_file = sys.argv[5]
    # Merged icefinder results
    icefinder_results_file = sys.argv[6]
    # Merged digis results
    digis_results_file = sys.argv[7]
    # Merged islandpath results
    islandpath_results_file = sys.argv[8]
    # Distance in bp between target MGE and AMR/VF
    distance = int(sys.argv[9])

    # read AMR and VF matrix to pandas dataframe
    amr_df = pd.read_csv(amr_matrix_file, sep='\t', header=0, index_col=0, keep_default_na=False)
    vf_df = pd.read_csv(vf_matrix_file, sep='\t', header=0, index_col=0, keep_default_na=False)

    # match AMR with MGE in the same assembly contig
    output_amr_mge_df = amr_df
    samples = list(amr_df.index)
    amrs = list(amr_df.columns)

    for i in range(len(samples)):
        for j in range(len(amrs)):
            sample = samples[i]
            pos = amr_df.iat[i, j]

            if not pos:
                output_amr_mge_df.iat[i, j] = pos
            else:
                # Prophage - Phigaro
                amr_prophage_list = match_mge(sample, pos, phigaro_results_file, distance)
                if amr_prophage_list:
                    output_amr_mge_df.iat[i, j] = "{}|(Prophage){}".format(pos, ','.join(amr_prophage_list))

                # ICE - ICEfinder
                amr_ice_list = match_mge(sample, pos, icefinder_results_file, distance)
                if amr_ice_list:
                    output_amr_mge_df.iat[i, j] += "|(ICE){}".format(','.join(amr_ice_list))

                # Plasmid - Plasmidfinder
                # amr_plasmidfinder_list = match_mge(sample, pos, plasmidfinder_results_file)
                # if amr_plasmidfinder_list:
                #     output_amr_mge_df.iat[i, j] += "|(Plasmidfinder){}".format(','.join(amr_plasmidfinder_list))

                # Plasmid - Platon
                # amr_platon_list = match_mge(sample, pos, platon_results_file)
                # if amr_platon_list:
                #     output_amr_mge_df.iat[i, j] += "|(Platon){}".format(','.join(amr_platon_list))

                # IS - digIS
                amr_is_list = match_mge(sample, pos, digis_results_file, distance)
                if amr_is_list:
                    output_amr_mge_df.iat[i, j] += "|(IS){}".format(','.join(amr_is_list))

                # GI - IslandPath-DIMOB
                amr_gi_list = match_mge(sample, pos, islandpath_results_file, distance)
                if amr_gi_list:
                    output_amr_mge_df.iat[i, j] += "|(GI){}".format(','.join(amr_gi_list))

    # Output matched AMR and MGE results
    output_amr_mge_df.to_csv(amr_mge_output, sep="\t", header=True, index=True)

    # match VF with MGE in the same assembly contig
    output_vf_mge_df = vf_df
    samples = list(vf_df.index)
    vfs = list(vf_df.columns)

    for i in range(len(samples)):
        for j in range(len(vfs)):
            sample = samples[i]
            pos = vf_df.iat[i, j]

            if not pos:
                output_vf_mge_df.iat[i, j] = pos
            else:
                # Prophage - Phigaro
                vf_prophage_list = match_mge(sample, pos, phigaro_results_file, distance)
                if vf_prophage_list:
                    output_vf_mge_df.iat[i, j] = "{}|(Prophage){}".format(pos, ','.join(vf_prophage_list))

                # ICE - ICEfinder
                vf_ice_list = match_mge(sample, pos, icefinder_results_file, distance)
                if vf_ice_list:
                    output_vf_mge_df.iat[i, j] += "|(ICE){}".format(','.join(vf_ice_list))

                # Plasmid - Plasmidfinder
                # vf_plasmidfinder_list = match_mge(sample, pos, plasmidfinder_results_file, distance)
                # if vf_plasmidfinder_list:
                #     output_vf_mge_df.iat[i, j] += "|(Plasmidfinder){}".format(','.join(vf_plasmidfinder_list))

                # Plasmid - Platon
                # vf_platon_list = match_mge(sample, pos, platon_results_file, distance)
                # if vf_platon_list:
                #     output_vf_mge_df.iat[i, j] += "|(Platon){}".format(','.join(vf_platon_list))

                # IS - digIS
                vf_is_list = match_mge(sample, pos, digis_results_file, distance)
                if vf_is_list:
                    output_vf_mge_df.iat[i, j] += "|(IS){}".format(','.join(vf_is_list))

                # GI - IslandPath-DIMOB
                vf_gi_list = match_mge(sample, pos, islandpath_results_file, distance)
                if vf_gi_list:
                    output_vf_mge_df.iat[i, j] += "|(GI){}".format(','.join(vf_gi_list))

    # Output matched VF and MGE results
    output_vf_mge_df.to_csv(vf_mge_output, sep="\t", header=True, index=True)
