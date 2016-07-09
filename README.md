# CERNIA
ceRNA Prediction Algorithm

############################################################################################

Copyright 2016 Rosalba Giugno

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

############################################################################################

CERNIA is a computational method that perform a classification of gene pairs in order to
predict whether they crosstalk through known, validated miRNA-target interactions.


################################
   **Usage**
################################

The software is written in R and require some packages (see section 'Dependencies') to be run.
It's structured in a main file and several subdirectories with scripts and data used by the code.

In order to run the example code, we created a main file called Cernia.R that loads the data
of validated pairs, the MREs and the DT-Hybrid's recommendations, the tissue-specific
gene expression and compute the vector of scores. Then, it's performed the classification
of all possible pairs using the 10 batches of 100 svms each, already trained.

The output produces a file that contain the scored pairs that have four tab-separated fields:

ceR1	ceR2	Score	Crosstalk (Yes/No)


################################
   **Dependencies**
################################

CERNIA require the following libraries installed in R:
 - parallel, for cpu parallel caomputation
 - e1071, for machine learning classification and analysis
 

################################
   **Contacts**
################################

For any help, or doubts, please contact: rosalba.giugno@univr.it