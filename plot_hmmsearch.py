import glob

import Bio.SearchIO.HmmerIO as hmmer
import numpy as np
import pandas as pd
import plotly.figure_factory as ff

gene_presence_list_all = []

lps_domains = ["LptD", "LptF_LptG", "Prenyltransf", "LpxD", "Hexapep", "Hexapep_2", "Acetyltransf_11", "LpxB",
               "Sigma70_r4", "Sigma70_r4_2", "Sigma70_ECF", "Epimerase", "CTP_transf_3", "Glyco_trans_1_4",
               "Glycos_transf_1", "LpxC", "CBS", "SIS", "Hydrolase_3", "Hydrolase", "LptD_N", "ABC_tran", "EptA_B_N",
               "Sulfatase",
               "PagL", "Glycos_transf_N", "Glycos_transf_1", "Kdo", "Glyco_transf_4",
               "Glyco_transf_9"]

fab_domains = ['ACCA', 'Acyl_transf_1', 'Thiolase_C', 'adh_short_C2', 'FabA','adh_short']
fab_names = ['ACCA [Acetyl-CoA carboxylase]', 'Acyl_transf_1 [malonyl-CoA:acyl carrier protein transacylase]', 'Thiolase_C [-ketoacyl-ACP synthas125.e - FabDBH]', 'adh_short_C2 [enoyl-ACP reductase - FabI]', 'FabA [FabZA]','adh_short [3-ketoacyl-acp reductase - FabG]']

# NAME  ACCA
# NAME  Acyl_transf_1
# NAME  FabA
# NAME  Ketoacyl-synt_C
# NAME  Thiolase_C
# NAME  adh_short
# NAME  adh_short_C2

lps_new = ['LpxB', 'LpxC', 'LpxD', 'Lip_A_acyltrans', 'Metallophos', 'PAP2', 'Hexapep', 'Acetyltransf_11', 'Hexapep_2',
           'Glycos_transf_1', 'Glycos_transf_N', 'LpxK']

pfam_domain_to_product_dict = {}

rbox_domains = ['CoA_trans', "3HCDH", "3HCDH_N", 'ECH_1',
                'ECH_2', "Acyl-CoA_dh_1", "Acyl-CoA_dh_2", "Acyl-CoA_dh_M", "Acyl-CoA_dh_N",
                'ETF_alpha', 'ETF', 'Thiolase_C', 'Thiolase_N', 'AcetylCoA_hyd_C', 'AcetylCoA_hydro', 'Rnf-Nqr']

rbox = ['CoA_trans (ydif)',"3HCDH (HBD) ", "3HCDH_N (HBD)",  'ECH_1 (CRT)',
        'ECH_2 (CRT)', "Acyl-CoA_dh_1 (ACDH)", "Acyl-CoA_dh_2 (ACDH)",
        "Acyl-CoA_dh_M (ACDH)",
        "Acyl-CoA_dh_N (ACDH)",
        'ETF_alpha (ETF-A)', 'ETF (ETF-B)', 'Thiolase_C (Thl C-terminal domain)',
        'Thiolase_N (Thl N-terminal domain)', 'Thioesterase_C', 'Thioesterase_hydro', 'RNF']


# RBOX PFAM MODEL
'''
NAME  3HCDH
NAME  3HCDH_N
NAME  Acyl-CoA_dh_1
NAME  Acyl-CoA_dh_2
NAME  Acyl-CoA_dh_M
NAME  Acyl-CoA_dh_N
NAME  ECH_1
NAME  ECH_2
NAME  ETF_alpha
NAME  ETF
NAME  Thiolase_C
NAME  Thiolase_N
NAME  AcetylCoA_hyd_C
NAME  AcetylCoA_hydro
NAME  CoA_trans
NAME  Rnf-Nqr
NAME  Thioesterase
'''

pfam_domain_to_product_dict["LptD"] = "LPS-assembly protein LptD"
pfam_domain_to_product_dict["LptF_LptG"] = "Lipopolysaccharide export system permease protein LptF/LptG"
pfam_domain_to_product_dict["Prenyltransf"] = "Ditrans,polycis-undecaprenyl-diphosphate synthase ((2E,6E)-farnesyl-diphosphate specific)"
pfam_domain_to_product_dict["LpxD"] = "UDP-3-O-(3-hydroxymyristoyl)glucosamine N-acyltransferase"
pfam_domain_to_product_dict["Hexapep"] = "Acyl-[acyl-carrier-protein]--UDP-N-acetylglucosamine O-acyltransferase / UDP-3-O-acylglucosamine N-acyltransferase 1"
pfam_domain_to_product_dict["Hexapep_2"] = "Acyl-[acyl-carrier-protein]--UDP-N-acetylglucosamine O-acyltransferase / UDP-3-O-acylglucosamine N-acyltransferase 1"
pfam_domain_to_product_dict["Acetyltransf_11"] = "Acyl-[acyl-carrier-protein]--UDP-N-acetylglucosamine O-acyltransferase"
pfam_domain_to_product_dict["LpxB"] = "Lipid-A-disaccharide synthase"
pfam_domain_to_product_dict["Sigma70_r2"] = "ECF RNA polymerase sigma-E factor"
pfam_domain_to_product_dict["Sigma70_r4"] = "ECF RNA polymerase sigma-E factor"
pfam_domain_to_product_dict["Sigma70_r4_2"] = "ECF RNA polymerase sigma-E factor"
pfam_domain_to_product_dict["Sigma70_ECF"] = "ECF RNA polymerase sigma-E factor"
pfam_domain_to_product_dict["Epimerase"] = "UDP-N-acetylglucosamine 4-epimerase"
pfam_domain_to_product_dict["Epimerase_2"] = "UDP-N-acetylglucosamine 4-epimerase"
pfam_domain_to_product_dict["CTP_transf_3"] = "3-deoxy-manno-octulosonate cytidylyltransferase"
pfam_domain_to_product_dict["Glyco_trans_1_4"] = "Glycogen(starch) synthase"
pfam_domain_to_product_dict["Glyco_transf_1"] = "Glycogen(starch) synthase"
pfam_domain_to_product_dict["LpxC"] = "UDP-3-O-acyl-N-acetylglucosamine deacetylase"
pfam_domain_to_product_dict["CBS"] = "Arabinose-5-phosphate isomerase"
pfam_domain_to_product_dict["SIS"] = "Arabinose-5-phosphate isomerase"
pfam_domain_to_product_dict["Hydrolase_3"] = "3-deoxy-manno-octulosonate-8-phosphatase"
pfam_domain_to_product_dict["Hydrolase"] = "3-deoxy-manno-octulosonate-8-phosphatase"
pfam_domain_to_product_dict["LptD_N"] = "Lipopolysaccharide export system protein LptA"
pfam_domain_to_product_dict["ABC_tran"] = "Hemin import ATP-binding protein HmuV"
pfam_domain_to_product_dict["EptA_B_N"] = "Phosphoethanolamine transferase CptA"
pfam_domain_to_product_dict["Sulfatase"] = "Phosphoethanolamine transferase CptA"
pfam_domain_to_product_dict["PagL"] = "Lipid A deacylase PagL"
pfam_domain_to_product_dict["Glycos_transf_N"] = "Lipid IV(A) 3-deoxy-D-manno-octulosonic acid transferase"
pfam_domain_to_product_dict["Glycos_transf_1"] = "Lipid IV(A) 3-deoxy-D-manno-octulosonic acid transferase"
pfam_domain_to_product_dict["Kdo"] = "Lipopolysaccharide core heptose(I) kinase RfaP"
pfam_domain_to_product_dict["Glycos_transf1"] = "Lipopolysaccharide core biosynthesis protein RfaG"
pfam_domain_to_product_dict["Glyco_transf_4"] = "Lipopolysaccharide core biosynthesis protein RfaG"
pfam_domain_to_product_dict["Glyco_transf_9"] = "Lipopolysaccharide heptosyltransferase 1"

pfam_description = [pfam_domain_to_product_dict[x] for x in lps_domains]

# test = [pfam_description[pfam_description]]
# print(pfam_description)
pfam_description = np.array(pfam_description)


# pfam_stacked
def get_query_results(hmmer_file):
    results = []
    file = open(hmmer_file, 'r')
    hmmer_parser = hmmer.hmmer3_domtab.Hmmer3TabParser(file)

    for query in hmmer_parser:
        print(query.id)
        # print(len(query))
        # print(query)
        results.append((query.id, len(query)))
        print(hmmer_file)

    # results = list(set(results))
    # print(results)
    # print(results)
    return results


def get_otu_name_from_filename(file_name: str):
    otu_name = file_name.split("/")[-1]
    otu_name = otu_name.split(".")[0]
    otu_name = otu_name.split("name")[1]
    return otu_name


# create dictionary with OTU as keys and list with query hits as values
gene_presence_dict = {}
gene_presence_dict = {}
gene_presence_dict_all = {}

def get_results_and_plot(pfam_domain_list, plot_path, plot_title, file_path):
    files = sorted(glob.glob(file_path))

    sum_genes = 0
    heatmap_text = []

    for file in files:
        if pfam_domain_list == "pfam_domains":
            otu = get_otu_name_from_filename(file)
        else:
            otu = file.split("/")[-1]
            otu = otu.split(".fasta")[0]

            # otu = otu.split("_")[1]
            # otu = otu.replace("-"," ")
        # fille dictionary with all hmmer search results in folder

        # get query results for all files and check if pfam id is present in query result - store result in list and then
        # dict for plotting

        gene_presence_list = []
        
        heatmap_annotation_one_file = []
        results = get_query_results(file)

        if len(results) == 0:

            continue

        # print("Results")
        # print(results)

        # check if otu has rbox
        for i, res in enumerate(results):
            print(f"i: {i} - res: {res}")
            # print(res)
            # print(f"i: {i} ------- len results {len(results)}")
            if res[0] != 'AcetylCoA_hyd_C' and res[0] != 'AcetylCoA_hydro':
                # print(res[1])
                # print(f"RESI:T SAIDJASPDJASPODJ {res}")
                if res[1] < 1:
                    continue




        # print(f"LENGTH OF HAS RBOX LIST: {len(has_rbox)}")

        # print(has_rbox)
        # print(list(gene_presence_dict.keys()))
        print(file_path)


        query_result_names = list(list(zip(*results))[0])
        query_result_counts = list(list(zip(*results))[1])
        for pdom in pfam_domain_list:
            for id, name in enumerate(query_result_names):
                count = 0

                if name == pdom:
                    count=count+query_result_counts[id]

                    if count >= 10:
                        count = 10



                    # gene_presence_list.append(1)

                    gene_presence_list.append(count)
                    gene_presence_list_all.append(count)


            if pdom not in query_result_names:
                gene_presence_list.append(0)
                gene_presence_list_all.append(0)


        heatmap_text.append(gene_presence_list)
        print(gene_presence_list)



        # print(heatmap_text)

        # for i,otu in enumerate(list(gene_presence_dict.keys())):
        #     for rbox_count in otu:
        #         print(rbox_count)

        sum_genes = sum(gene_presence_list)
        # print(sum_genes)
        # skip file if gene_presence list is empty = no query hits for for otu
        if sum_genes < 5:
            continue

        gene_presence_dict[otu] = gene_presence_list
        gene_presence_dict_all[otu] = gene_presence_list_all

    # print(f"Total sum of gene counts: {sum_genes}")
    df = pd.DataFrame.from_dict(gene_presence_dict)
    df = df.to_numpy()

    df_all = pd.DataFrame.from_dict(gene_presence_dict_all)
    df_all = df_all.to_numpy()

    # rbox = ["3HCDH (HBD) ", "3HCDH_N (HBD)", "Acyl-CoA_dh_1 (ACDH)", "Acyl-CoA_dh_2 (ACDH)", "Acyl-CoA_dh_M (ACDH)", "Acyl-CoA_dh_N (ACDH)", 'ECH_1 (CRT)',
    #             'ECH_2 (CRT)',
    #             'ETF_alpha (ETF-A)', 'ETF (ETF-B)', 'Thiolase_C (Thl C-terminal domain)', 'Thiolase_N (Thl N-terminal domain)', 'AcetylCoA_hyd_C (YDIF C-terminal domain)', 'AcetylCoA_hydro (YDIF N-terminal domain)', 'CoA_trans (ydif)']


    lps = ["LptD (LPS-assembly protein)", "LptF_LptG (Lipopolysaccharide export system permease protein LptF + LptG)",
           "Prenyltransf (uppS)", "LpxD", "Hexapep (lpxA)", "Hexapep_2 (lpxA)", "Acetyltransf_11 (lpxA)", "LpxB",
           "Sigma70_r4 (ECF RNA polymerase sigma-E factor", "Sigma70_r4_2 (ECF RNA polymerase sigma-E factor",
           "Sigma70_ECF (ECF RNA polymerase sigma-E factor)",
           "Epimerase (wbpP)", "CTP_transf_3 (kdsB)", "Glyco_trans_1_4 (Glycosyl transferase)",
           "Glycos_transf_1 (Glycosyl transferase)", "LpxC", "CBS (kdsD|kpsF)", "SIS (kdsD|kpsF)",
           "Hydrolase_3 (kdsC - 3-deoxy-manno-octulosonate-8-phosphatase activity) ",
           "kdsC - Hydrolase (3-deoxy-manno-octulosonate-8-phosphatase activity)",
           "LptD_N (Lipopolysaccharide export system protein LptA)", "ABC_tran (Hemin import ATP-binding protein HmuV)",
           "EptA_B_N (eptC - Phosphoethanolamine transferase CptA)",
           "Sulfatase (eptC - Phosphoethanolamine transferase CptA)",
           "PagL (Lipid A deacylase PagL)",
           "Glycos_transf_N (Lipid IV(A) 3-deoxy-D-manno-octulosonic acid transferase)",
           "Glycos_transf_1 (Lipid IV(A) 3-deoxy-D-manno-octulosonic acid transferase)",
           "Kdo (Lipopolysaccharide core heptose(I) kinase RfaP)",
           "Glyco_transf_4 (Lipopolysaccharide core biosynthesis protein RfaG)",
           "Glyco_transf_9 (Lipopolysaccharide heptosyltransferase 1)"]

    lps_names_new = ['LpxB', 'LpxC', 'LpxD', 'Lip_A_acyltrans (LpxL + LpxM)', 'Metallophos (LpxH)',
                     'PAP2 (lpexe + lpxf)', 'Hexapep (LpxA 1)', 'Acetyltransf_11 (LpxA 2)', 'Hexapep_2 (LpxA 3)',
                     'Glycos_transf_1 (kdta 2)', 'Glycos_transf_N (kdta 1)', 'LpxK']
    methanogen_names_new = ['ADH_N - Alcohol dehydrogenase', 'ADH_N_assoc - Alcohol dehydrogenase',
                            'Aldedh - aldehyde dehydrogenase']


    b_oxidation_names = ['Acyl-CoA_dh_1 [acyl-coa dehydrogenase]', 'Acyl-CoA_dh_M [acyl-coa dehydrogenase]', 'Acyl-CoA_dh_N [acyl-coa dehydrogenase]', 'ECH_1 [enyol coa hydratase]', 'Thiolase_C [beta-ketothiolase]', 'Thiolase_N [beta-ketothiolase]',
                         'adh_short [3-hydroxyacyl-coa dehydrogenase]']

    if pfam_domain_list == lps_domains:
        pfam_domain_list = lps
    if pfam_domain_list == rbox_domains:
        pfam_domain_list = rbox
    if pfam_domain_list == lps_new:
        pfam_domain_list = lps_names_new
    # if pfam_domain_list == methanogen_domains:
    #     pfam_domain_list = methanogen_names_new
    if pfam_domain_list == fab_domains:
        pfam_domain_list = fab_names
    # if pfam_domain_list == b_oxidation_domains:
    #     pfam_domain_list = b_oxidation_names

    fig = ff.create_annotated_heatmap(z=df, y=pfam_domain_list, x=list(gene_presence_dict.keys()))
    print(gene_presence_dict.keys())
    heatmap_all = ff.create_annotated_heatmap(z=df, y=pfam_domain_list, x=list(gene_presence_dict.keys()))

    # fig.add_trace(go.Heatmap(z=df, y=pfam_domain_list, x=list(gene_presence_dict.keys()),zmin=0,colorscale="blackbody"))
    fig.update_layout(title=plot_title, xaxis=dict(side='bottom'), width=2000,
                      height=1300)

    heatmap_all.update_layout(title=plot_title, xaxis=dict(side='bottom'), width=2000,
                      height=1300)


    # print(list(gene_presence_dict.keys()))
    # fig is the heatmap for each file... heatmap_all is the heatmap for all files... the last plot that is built for
    # heatmap all is the final file. the other ones can be ignored and are the plots that show the genes from n=1..9
    # where n is the index of the reactor
    # plotly.offline.plot(fig, filename=plot_path, auto_open=False)
    # fig.update_xaxes(side="top")
    # fig.show()
    heatmap_all.update_xaxes(side="bottom")
    heatmap_all.show()






# get_results_and_plot(fab_domains, "/Users/timolucas/Documents/spirito/assembly/r1t1/fab_r1.html",
#                      "FAB domains r1",
#                      '/Users/timolucas/Documents/spirito/assembly/r1t1/fab/*.txt')
#

# get_results_and_plot(fab_domains, "/Users/timolucas/Documents/spirito/assembly/r2t1/fab_r2.html",
#                      "FAB domains r2",
#                      '/Users/timolucas/Documents/spirito/assembly/r2t1/fab/*.txt')
# #
# get_results_and_plot(fab_domains, "/Users/timolucas/Documents/spirito/assembly/r3t1/fab_r3.html",
#                      "FAB domains r3",
#                      '/Users/timolucas/Documents/spirito/assembly/r3t1/fab/*.txt')
# #
# #


# get_results_and_plot(rbox_domains, "/Users/timolucas/Documents/spirito/proteomics/rbox.html",
#                      "RBOX domains proteomic reads",
#                      '/Users/timolucas/Documents/spirito/proteomics/final/rbox/*.txt')


# get_results_and_plot(fab_domains, "/Users/timolucas/Documents/spirito/proteomics/fab.html",
#                      "FAB domains proteomic reads sept 2021",
#                      '/Users/timolucas/Documents/spirito/proteomics/final/fab/*.txt')
#

# get_results_and_plot(fab_domains, "/Users/timolucas/Documents/spirito/assembly/r3t1/fab_r3.html",
#                      "FAB domains r3",
#                      '/Users/timolucas/Documents/spirito/assembly/r3t1/fab/*.txt')



# get_results_and_plot(b_oxidation_domains, "/Users/timolucas/Documents/jeon/r1_august/b_oxidation_r1_binary.html",
#                      "B_oxidation proteins in reactor 1",
#                      '/Users/timolucas/Documents/jeon/r1_august/prokka/b_oxidation/*.txt')


# get_results_and_plot(rbox_domains, "/Users/timolucas/Documents/spirito/proteomics/rbox_ncbi.html",
#                      "RBOX domains proteomic reads NCBI names",
#                      '/Users/timolucas/Documents/spirito/proteomics/ncbi/rbox/*.txt')

# get_results_and_plot(fab_domains, "/Users/timolucas/Documents/spirito/proteomics/fab_ncbi.html",
#                      "FAB domains proteomic reads NCBI names",
#                      '/Users/timolucas/Documents/spirito/proteomics/ncbi/fab/*.txt')

fab_paths=["/Users/timolucas/Documents/spirito/hmmsearch/r1t1/fab","/Users/timolucas/Documents/spirito/hmmsearch/r1t2/fab",
           "/Users/timolucas/Documents/spirito/hmmsearch/r1t3/fab","/Users/timolucas/Documents/spirito/hmmsearch/r2t1/fab",
           "/Users/timolucas/Documents/spirito/hmmsearch/r2t2/fab","/Users/timolucas/Documents/spirito/hmmsearch/r2t3/fab",
           "/Users/timolucas/Documents/spirito/hmmsearch/r3t1/fab","/Users/timolucas/Documents/spirito/hmmsearch/r3t2/fab",
           "/Users/timolucas/Documents/spirito/hmmsearch/r3t3/fab"]


rbox_paths=["/Users/timolucas/Documents/spirito/hmmsearch/r1t1/rbox","/Users/timolucas/Documents/spirito/hmmsearch/r1t2/rbox",
           "/Users/timolucas/Documents/spirito/hmmsearch/r1t3/rbox","/Users/timolucas/Documents/spirito/hmmsearch/r2t1/rbox",
           "/Users/timolucas/Documents/spirito/hmmsearch/r2t2/rbox","/Users/timolucas/Documents/spirito/hmmsearch/r2t3/rbox",
           "/Users/timolucas/Documents/spirito/hmmsearch/r3t1/rbox","/Users/timolucas/Documents/spirito/hmmsearch/r3t2/rbox",
           "/Users/timolucas/Documents/spirito/hmmsearch/r3t3/rbox"]

#metagenomics plots fab

for path in fab_paths:
    sample=path.split("/")[-2]
    get_results_and_plot(fab_domains, f"{path}/{sample}_fab.html",
                         f"FAB domains in metagenome sample {sample}",
                         f"{path}/*.txt")

# for path in rbox_paths:
    # sample=path.split("/")[-2]
    # get_results_and_plot(rbox_domains, f"{path}/{sample}_rbox.html",
    #                      f"RBOX domains in metagenome after adding counts of sample {sample}",
    #                      f"{path}/*.txt")




#
# ['aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance',
#  'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg',
#  'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl',
#  'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric',
#  'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys',
#  'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet',
#  'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges',
#  'orrd', 'oryel', 'peach', 'phase', 'picnic', 'pinkyl', 'piyg',
#  'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn', 'puor',
#  'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu', 'rdgy',
#  'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar', 'spectral',
#  'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn', 'tealrose',
#  'tempo', 'temps', 'thermal', 'tropic', 'turbid', 'twilight',
#  'viridis', 'ylgn', 'ylgnbu', 'ylorbr', 'ylorrd'].


# has RBOX R1
# [False, False, False, False, False, True, False, False, False, True, True, False, False, False, True, True, False, False, True, True, False, False, False, False, True, False, True, False, False, True, True, False, False, False, False, False, False, False, True, False, True, False, True, True, False, True, False, True, False, True, False, True, False, False, False, False, False, False, False, False, False, False, True, True, False, False, False, False, False, False, True, False, False, False, False, False, False, True, False, False, False, False, False, False, False, True, False, False, True, True, False, False, False, False, True, False, False, False, False, False, False, False, True, False, True, False, True, False, False, False, False, False, False, True, True, True, False, True, True, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, True, False, False, True, False, False, True, False, False, True, True, False, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, True, True, False, False, True, True, True, True, False, False, True, False, False, False, False, False, True, False, False, False, True, False, False, False, True, True, False, True, False, True, False, True, False, False, True, False, False, True, True, False, False, True, False, False, False, False, False, True, False, True, True, False, False, True, False, True, False, True, False, False, False, False, True, False, True, False, False, False, True, True, False, False, False, False, True, True, False, True, False, False, False, False, False, True]
