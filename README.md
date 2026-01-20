# Evolutionary bioenergetics of sporulation

**Last updated:** 19 November 2025  
**Authors:** Canan Karakoç,  William R. Shoemaker, Jay T. Lennon

This repository contains the complete analysis pipeline of the manuscript "Evolutionary bioenergetics of sporulation". 
The project integrates genomic, proteomic, transcriptomic, and by-the-number datasets to quantify ATP investment across the spore lifecycle, 
compares with other cellular investments.
It also incorporates these estimates into mechanistic models to show how metabolic constraints shape sporulation efficiency, and how this developmental program becomes visible to selection. 

---

## Project Structure

- bioaccounting/:   Bioenergetic calculations using by-the-numbers, spore formation and revival data.
- efficiency/:      Data and plots related to sporulation efficiency across studies. 
                  Models for population dynamics in batch and chemostat conditions and sporulation efficiency.
- evolution/:       Data and figures of sporulation COGs. Model outputs and figures related to evolutionary dynamics.


---

## Getting Started

- All by-the-numbers calculations in the manuscript are embedded in the  R code with their references, 
such as membrane area and membrane costs calculations, genome size, nucleotide costs, and replication costs.
- Key datasets such as amino acid costs, spore-formation, and revival genes/proteins are listed below. 
- The project includes models and simulations written in Python, and the codes are provided in their respective folder. 

To reproduce the analysis:

1. Run `code\cost_of_spore.R` from the top for bioaccounting. Ensure the working directory is set properly.
2. Run Python codes to reproduce data for population dynamics, spore efficiency, and evolutionary analysis. 
Codes are stored in "efficiency" and "evolution" sections under "model" folders. 
3. Run `efficiency_plots.R` and `COGs_figures.R` to reproduce figures of empirical and model data of 
population dynamics, efficiency, and evolutionary outcomes. These codes are stored under respective sections and "empirical" folders. 

---

## Key Datasets Used for Bioaccounting 

| Dataset Description | Source |
|---------------------|--------|
| Gene & protein annotation | [SubtiWiki](https://subtiwiki.uni-goettingen.de/) |
| Expression data | [SporeWeb](https://sporeweb.molgenrug.nl/) |
| Protein abundance | [PAX-db](https://pax-db.org/) |
| Newly synthesized proteins during germination | [Swarge et al. 2020](https://doi.org/10.1128/mSphere.00463-20) |
| Protein sequences | [UniProt: B. subtilis 168](https://www.uniprot.org/taxonomy/224308) |
| Amino acid/nucleotide ATP costs | [Mahmoudabadi et al. 2017](https://doi.org/10.1073/pnas.1701670114) |

---

## Cost definitions used throughout the analysis

| Cost Type | Description |
|-----------|-------------|
| **Opportunity Cost** | ATP-equivalent energy required to synthesize building blocks. |
| **Direct Cost**      | ATPs used in polymerization and structural assembly. |
| **Total Cost**       | Sum of opportunity + direct costs. |

---

## Empirical data for COGs and efficiency 

Spore-related COGs are derived from [Galperin et al. 2022](https://journals.asm.org/doi/10.1128/jb.00079-22). 
See manuscript, supplementary information for empirical efficiency data derived from various sources. 

---

## Main Outputs

- **Figure 1**: Cumulative ATP costs of sporulation over time. (bioaccounting/)
- **Figure 2**: Energetic costs of germination and outgrowth. (bioaccounting/)
- **Figure 3**: Comparison with alternative stress responses and cellular processes. (bioaccounting/)
- **Figure 4**: Empirical and theoretical models of sporulation efficiency. (efficiency/)
- **Figure 5**: COG distributions and evolutionary constraints. (evolution/)

---


##  References

See embedded DOIs in script comments. Full reference list available in the manuscript.

---

##  Reproducibility Notes

- All scripts are annotated and reproducible within the RStudio project.
- Final manuscript figures formatted in Adobe Illustrator.  

---

## Contact

Jay T. Lennon — [lennonj@iu.edu](mailto:lennonj@iu.edu)
