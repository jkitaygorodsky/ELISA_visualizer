# ELISA_visualizer
This R Shiny tool has been developed for visualizing antibody levels associated with COVID-19, compatible with the summary output files provided by the Serology team at the Network Biology Collaborative Centre (NBCC) at Sinai Health.
Note, that this app is intended for use by clinical/research study staff ONLY, not for patient use. This app is not responsible for storage or sharing of any confidential patient information.

## Getting Started
Prerequisites
To run the tool, you need R version 3.6.3 and the following libraries:

* shiny
* shinymanager
* ggplot2
* tidyverse
* grid
* gtable
* gridExtra
* gridtext

Clone the repository:
```
git clone https://github.com/jkitaygorodsky/ELISA_visualizer
```

## Input
Prior to running the tool, ensure that the column names of your file are compatible with the app, and that Timepoint information is filled in.

* Spike antibody levels: SmT1 IgG [dilution], SmT1 BAU/ml aggregated
* RBD antibody levels: RBD IgG [dilution], RBD BAU/ml aggregated
* NP antibody levels: NP IgG [dilution], NP BAU/ml aggregated

## Running the Tool

To run the tool, open the app.R file, and run the app (e.g. in RStudio).

The tool will prompt you to load a .CSV file - this is the summary output file of the ELISA test results.

If a participant ID is entered, that participant's data points will be highlighted with asterisks on the graph, and the associated table will be populated with that participant's antibody levels across all measured timepoints.

Individual timepoints and dilutions may be filtered.

Data may be log2-transformed with the provided checkbox.

## Output
The graph, table, and datapoint annotations can be downloaded as a pdf.
