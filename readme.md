
# CellOntologyMapper:  Consensus mapping of cell type annotation

Single-cell RNA sequencing has revolutionized cellular biology, with atlases now encompassing over 100 million cells. However, researchers employ vastly different naming conventions when annotating cell types, creating a fragmented landscape that severely impedes data integration and comparative analysis. Here, we present CellOntologyMapper, an automated framework that standardizes cell type annotations by intelligently mapping user-defined names to established Cell Ontology and Cell Taxonomy identifiers. Our approach leverages advanced natural language processing, including sentence transformers and large language models, to interpret diverse naming conventions and resolve them to standardized ontological terms. The system handles complex challenges including abbreviated cell names, synonym resolution, and context-dependent interpretation. We built a comprehensive query system based on 19,381 cell type entries from established Cell Ontology and Cell Taxonomy databases, which organize into 24 biologically coherent clusters. Systematic validation across datasets spanning different scales—from 17-cell-type lung studies to 91-cell-type immune atlases and complex developmental systems—consistently demonstrated robust performance and high accuracy. CellOntologyMapper successfully resolved annotation challenges across conventional tissue studies, rare cell populations, and developmentally intricate datasets rich in abbreviated nomenclature. By providing automated, scalable annotation harmonization, our framework enables researchers to leverage existing single-cell datasets while ensuring compatibility with future atlas efforts. When incorporated into publications, CellOntologyMapper enhances the credibility and reference value of cell type annotations, representing a crucial step toward an integrated single-cell genomics ecosystem.

## Installation

```shell
pip install git+https://github.com/Starlitnightly/CellOntologyMapper.git
```

or you can direct use CellOntologyMapper with [OmicVerse](https://github.com/Starlitnightly/omicverse/)

The detailed tutorial could be found in https://omicverse.readthedocs.io/en/latest/Tutorials-single/t_cellmatch/

## Contact

- Zehua Zeng ([starlitnightly@gmail.com](mailto:starlitnightly@gmail.com) or [zehuazeng@xs.ustb.edu.cn](mailto:zehuazeng@xs.ustb.edu.cn))

## Other

If you would like to sponsor the development of our project, you can go to the afdian website (https://ifdian.net/a/starlitnightly) and sponsor us.


Copyright © 2025 [112 Lab](https://112lab.asia/). <br />
This project is [GPL3.0](./LICENSE) licensed.
