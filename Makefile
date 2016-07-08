USER = $(shell whoami)
MAKE = /usr/bin/make
DATA_DIR = data/processed
SRC = src
PLT_DIR = plotting_scripts
TEX_DIR = /Users/$(USER)/Dropbox/g1_leaf_ecosystem_paper/figures
FIG_DIR = /Users/$(USER)/Dropbox/g1_leaf_ecosystem_paper/figures/figs


all: data paper
data: $(DATA_DIR)/g1_leaf_gas_exchange.csv \
	  $(DATA_DIR)/g1_fluxnet.csv \
	  $(DATA_DIR)/g1_fluxnet_screened.csv \
	  $(DATA_DIR)/g1_fluxnet_PM.csv \
	  $(DATA_DIR)/g1_fluxnet_screened_PM.csv \
	  $(DATA_DIR)/g1_isotope_screened.csv \
	  $(DATA_DIR)/g1_fluxnet_screened_matching_sites.csv \
	  $(DATA_DIR)/g1_signifletters_withinPFT.csv \
	  $(DATA_DIR)/g1_siglet_withinPFT_coupled.csv
paper: $(TEX_DIR)/figures.pdf

##
# Make data files
##
$(DATA_DIR)/g1_leaf_gas_exchange.csv:	$(SRC)/estimate_g1_from_leaf_gas_exchange.py
	python $<

$(DATA_DIR)/g1_fluxnet_screened.csv:	$(SRC)/screen_fluxnet_g1_fits.py \
										$(DATA_DIR)/g1_fluxnet.csv
	python $<

$(DATA_DIR)/g1_fluxnet_PM.csv:	$(SRC)/estimate_g1_from_fluxnet_data_PM_MPI.py
	python $<

$(DATA_DIR)/g1_fluxnet_screened_PM.csv:	$(SRC)/screen_fluxnet_g1_fits.py \
										$(DATA_DIR)/g1_fluxnet_PM.csv
	python $<

$(DATA_DIR)/g1_isotope_screened.csv:	$(SRC)/estimate_g1_from_leaf_isotope.py
	python $<

$(DATA_DIR)/g1_signifletters_withinPFT.csv:	stats/stats_g1_withinPFT.R \
											$(DATA_DIR)/g1_fluxnet_screened_PM.csv \
											$(DATA_DIR)/g1_isotope_screened.csv \
											$(DATA_DIR)/g1_leaf_gas_exchange.csv
	R CMD BATCH stats/stats_g1_withinPFT.R /dev/null

$(DATA_DIR)/g1_fluxnet_screened_matching_sites.csv:	$(SRC)/screen_fully_coupled_fluxnet_to_match_inverted_sites.py \
													$(DATA_DIR)/g1_fluxnet_PM.csv \
													$(DATA_DIR)/g1_fluxnet_screened_PM.csv
	python $<

$(DATA_DIR)/g1_siglet_withinPFT_coupled.csv:	stats/stats_g1_withinPFT_fullycoupled.R \
												$(DATA_DIR)/g1_fluxnet_screened_PM.csv \
												$(DATA_DIR)/g1_isotope_screened.csv \
												$(DATA_DIR)/g1_leaf_gas_exchange.csv \
												$(SRC)/screen_fully_coupled_fluxnet_to_match_inverted_sites.py
	R CMD BATCH stats/stats_g1_withinPFT_fullycoupled.R /dev/null


##
# Figures - set these up so that if any data changes, they are *ALL* remade
##
$(FIG_DIR)/g1_gas_exchange_isotope_fluxnet_boxplot.pdf:	$(PLT_DIR)/plot_g1_boxplot.py \
														$(DATA_DIR)/*.csv
	python $<



$(FIG_DIR)/g1_vs_latitude.pdf:	$(PLT_DIR)/plot_g1_vs_latitude.py \
								$(DATA_DIR)/*.csv
	python $<

$(FIG_DIR)/g1_matching_leaf_and_canopy.pdf:	$(PLT_DIR)/plot_overlapping_leaf_flux.py \
											$(DATA_DIR)/*.csv
	python $<

$(FIG_DIR)/g1_vs_c4frac.pdf:	$(PLT_DIR)/plot_g1_vs_c4frac.py \
								$(DATA_DIR)/*.csv
	python $<

$(FIG_DIR)/g1_vs_lai.pdf:	$(PLT_DIR)/plot_g1_vs_lai.py \
							$(DATA_DIR)/*.csv
	python $<

$(FIG_DIR)/g1_gas_exchange_isotope_fluxnet_boxplot_inverted_PM.pdf:	$(PLT_DIR)/plot_g1_boxplot_inverted_PM.py \
																	$(DATA_DIR)/*.csv
	python $<

$(FIG_DIR)/location_map.png:	$(PLT_DIR)/make_map_of_site_locations.py \
								$(DATA_DIR)/*.csv
	python $<

# figure.pdf
##
$(TEX_DIR)/figures.pdf:	$(TEX_DIR)/figures.tex \
						$(FIG_DIR)/g1_gas_exchange_isotope_fluxnet_boxplot.pdf \
						$(FIG_DIR)/g1_vs_latitude.pdf \
						$(FIG_DIR)/g1_matching_leaf_and_canopy.pdf \
						$(FIG_DIR)/g1_vs_c4frac.pdf \
						$(FIG_DIR)/g1_vs_lai.pdf \
						$(FIG_DIR)/g1_gas_exchange_isotope_fluxnet_boxplot_inverted_PM.pdf \
						$(FIG_DIR)/location_map.png
	cd $(TEX_DIR) && $(MAKE) clean && $(MAKE)

cleanfigs:
	rm -f $(FIG_DIR)/*.pdf

.PHONY : clean

clean:
	rm -f $(DATA_DIR)/*.csv $(FIG_DIR)/*.pdf
