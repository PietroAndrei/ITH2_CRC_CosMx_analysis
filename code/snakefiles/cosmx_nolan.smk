configfile: "config.yaml"

rule all_cosmx_nolan:
	input:
		"../data/CosMx_Nolan/CosMx_data_Fig6I_CRC.rds",
		"intermediate_outputs/CosMx_Nolan/InSituType/unsup/Cosmx_Nolan_Unsup_Clust_res.rds",
		"intermediate_outputs/CosMx_Nolan/InSituType/unsup/Cosmx_Nolan_Unsup_FinalClust_res.rds",
		"../plots/FigS4D_FDR_barplot_TAMT_EpiT_APP_IFN_genes.pdf"

rule prepare_cosmx_Nolan_validation_data:
	output:
		cosmx_clean = "../data/CosMx_Nolan/CosMx_data_Fig6I_CRC.rds"
	params:
		data_path = "../data/CosMx_Nolan"
	conda:
		"cosmx"
	resources:
		mem_mb = 12000,
		mem_mb2 = 12000,
		runtime_min = 180,
		disk_mb = 1000
	script:
		config["scripts"]+"prepare_cosmx_Nolan_validation_data.R"

#Perform clustering analysis with InSituType on all CRC single cells
rule cosmx_nolan_insitutype_clustering:
	input:
		seurat = "../data/CosMx_Nolan/CosMx_data_Fig6I_CRC.rds"
	output:
		plots = "intermediate_outputs/CosMx_Nolan/InSituType/plots/CosMx_Nolan_InSituClust_UMAP_plus_FOVplot.pdf",
		markers = "intermediate_outputs/CosMx_Nolan/InSituType/unsup/DE_res/Cosmx_Nolan_MainClusters_markers.rds",
		cell_anno = "intermediate_outputs/CosMx_Nolan/InSituType/unsup/CellAnno/Cosmx_Nolan_MainClusters_cell_annotation.rds",
		InSitu_unsup_res = "intermediate_outputs/CosMx_Nolan/InSituType/unsup/Cosmx_Nolan_Unsup_Clust_res.rds"
	params:
		data_path = "../data/CosMx_Nolan",
		utils = "src/cosmx_utils.R"
	conda:
		"cosmx"
	resources:
		mem_mb = 15000,
		mem_mb2 = 15000,
		runtime_min = 360,
		disk_mb = 3000
	script:
		config["scripts"]+"cosmx_nolan_insitutype_clustering.R"

rule refine_cosmx_nolan_clusters:
	input:
		seurat = "../data/CosMx_Nolan/CosMx_data_Fig6I_CRC.rds",
		markers = "intermediate_outputs/CosMx_Nolan/InSituType/unsup/DE_res/Cosmx_Nolan_MainClusters_markers.rds",
		anno = "intermediate_outputs/CosMx_Nolan/InSituType/unsup/CellAnno/Cosmx_Nolan_MainClusters_cell_annotation.rds",
		unsup = "intermediate_outputs/CosMx_Nolan/InSituType/unsup/Cosmx_Nolan_Unsup_Clust_res.rds"
	output:
		new_anno = "intermediate_outputs/CosMx_Nolan/InSituType/unsup/CellAnno/Cosmx_Nolan_RefinedClusters_cell_annotation.rds",
		new_markers = "intermediate_outputs/CosMx_Nolan/InSituType/unsup/DE_res/Cosmx_Nolan_RefinedClusters_markers.rds",
		final_unsup = "intermediate_outputs/CosMx_Nolan/InSituType/unsup/Cosmx_Nolan_Unsup_FinalClust_res.rds",
		plots = "intermediate_outputs/CosMx_Nolan/InSituType/plots/CosMx_Nolan_InSitu_FinalClust_FOVplot.pdf",
		meta_sf = "../data/CosMx_Nolan/CosMx_data_Fig6I_CRC_metadata.rds"
	params:
		data_path = "../data/CosMx_Nolan",
		utils = "src/cosmx_utils.R"
	conda:
		"cosmx"
	resources:
		mem_mb = 15000,
		mem_mb2 = 15000,
		runtime_min = 360,
		disk_mb = 3000
	script:
		config["scripts"]+"refine_cosmx_nolan_clusters.R"

rule compare_epi_tam_markers_tcell_distance:
	input:
		seurat = "../data/CosMx_Nolan/CosMx_data_Fig6I_CRC.rds",
		meta = "../data/CosMx_Nolan/CosMx_data_Fig6I_CRC_metadata.rds"
	output:
		wilcox_res_table = "intermediate_outputs/CosMx_Nolan/CosMx_Nolan_CRC_TAM_Epi_wilcoxon_results_AllMarkers.rds",
		fdr_plot = "../plots/FigS4D_FDR_barplot_TAMT_EpiT_APP_IFN_genes.pdf"
	params:
		data_path = "../data/CosMx_Nolan",
		utils = "src/cosmx_utils.R"
	conda:
		"cosmx"
	resources:
		mem_mb = 10000,
		mem_mb2 = 10000,
		runtime_min = 60,
		disk_mb = 1000
	script:
		config["scripts"]+"compare_epi_tam_markers_tcell_distance_cosmx_nolan.R"
