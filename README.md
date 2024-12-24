# cpu-based-ngs-workflow
A computational pipeline designed for automated processing of Next-Generation Sequencing (NGS) data using CPU resources. 

## functionalities 

- dynamic directory building 
- passing userID and jobID through the pipeline 
- logging each step and error handling 
- skipping each step if it previously ran (for killed runs)
- flag to skip the QC step
- conditionally start the pipeline based on the kind of inputs
- conditionally alternative variant filtering and parsing



