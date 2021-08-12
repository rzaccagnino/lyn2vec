# lyn2vec
A feature embedding method for biological sequences

## Pre-requisites

Most of the methods are implemented in Python using `scikit-learn` (and utility
functions from `numpy` and `pandas`).
For ease of reproduction, the exact versions of the packages can be installed
through the **conda package manager** by installing
[`miniconda`](https://docs.conda.io/en/latest/miniconda.html) and replicating
the environment using the file `spec-file.txt` in this repository with the
command:

```
conda create --name lyndonenv --file spec-file.txt
```

## Data

Due to size constraints, data is available at [this
URL](https://drive.google.com/drive/folders/1_E-wKUA6PNSMqIa2jBGyFMiak4ilDsK3).
The dataset is composed of the following two files:

- `transcripts_genes.fa.gz`, containing the sequences of the transcripts used
  for training the model and that should be placed (and `gunzip`-ed) in the
  `training` directory.
- `sample_10M_genes.fastq.gz`, containing a sample of 10M 100-bp long reads
  simulated from Human chromosomes 1, 17, and 21 using [Flux
  Simulator](https://dx.doi.org/10.1093/nar/gks666). This file should be placed
  (and `gunzip`-ed) in the `testing` directory.


## Reproduction

To carry out an experiment, follows the steps:

- `Compute fingerprints` 
    ```
    - `method`   : experiment_fingerprint_1f_np_step in SCRIPT fingerprint.py

    - `cmd line` : python fingerprint.py --type 1f_np --path training/ 
                 --fasta transcripts_genes.fa --type_factorization ICFL_COMB 
                 -n 8

    - `return`   : for type_factorization, computes a "fingerprint" file containing 
                   a row for each read, with the format "IDGENE FINGERPRINT", where 
                   "FINGERPRINT" is the fingerprint of the read
    ```
    
    `N.B.`
    ```
    - `--fact no_create` : do not create a file containing the factors corresponding 
                      to the fingerprint fingerprint (default: --fact create)
    - `--shift no_shift` : do not generate the shifts of lengths 100 (default: --shift shift)
    - `--filter no_list` : do not consider only the reads for the genes contained in 
                      the file list_experiment.txt (default: --filter no_list)
    ```

- `Build datasets`
    ```
    - `method`   : experiment_dataset_step in SCRIPT training_mp.py

    - `cmd line` : python training.py --step dataset --path training/ 
                   --type_factorization ICFL_COMB  --k_value 5 n 8

    - `return`   : for type_factorization it uses the corresponding 
                   "fingerprint" file to generate a dataset for each value of k. 
                   Such a dataset will be splitted in 2 pickle files: 
                   dataset_X_factorization which contain the samples, and 
                   dataset_y_factorization which contains the corresponding labels
    ```
    
    
- `Train K-fingers classifiers`
    ```
    - `method`   : experiment_training_step in SCRIPT training.py

    - `cmd line` : python training.py --step train --path training/ --k_value 5
                 --type_factorization ICFL_COMB  --model RF -n 4

    - `return`   : for each trained classifier save a PICKLE file 
                 (ex. RF_ICFL_COMB_K5.pickle) and the report CSV containing 
                 the metrics for the performance in training 
                 (ex. RF_kfinger_clsf_report_ICFL_COMB_K5.csv)
    ```

- `Reads classification`

    `pre-settings:`
    ```
    - A k-finger trained classifier (ex. RF_ICFL_COMB_K5.pickle)
    - The dataset for the k-finger trained classifier chosen 
       (ex. dataset_X_ICFL_COMB_K5.pickle, dataset_y_ICFL_COMB_K5.pickle)
    - The fingerprint and fact_fingerprint corresponding to the type of 
       factorization for which the chosen classifier was trained 
       (ex. fingerprint_ICFL_COMB.txt e fact_fingerprint_ICFL_COMB.txt)
    ```
        
    `RF Fingerprint classifier:` 
        
        `Training:`
        ```
          - `method`    : training_train_RF_fingerprint_step in SCRIPT training.py

          - `cmd line`  : python training.py --step train_RF_fingerprint 
                        --path training/ --model RF_FINGERPRINT
                        --type_factorization ICFL_COMB -n 4

          - `return`    : save the PICKLE RF FINGERPRINT trained 
                        (ex. RF_fingerprint_classifier_ICFL_COMB.pickle)
                        and the corresponding CSV report 
                        (RF_fingerprint_clsf_report_ICFL_COMB.csv")
        ```    
        
        `Testing reads:`
        ``` 
          - `method`    : testing_reads_RF_fingerprint_step in SCRIPT testing.py

          - `cmd line`  : python testing.py --step test_RF_fingerprint --path testing/ 
                          --rf_fingerprint_model RF_fingerprint_classifier_ICFL_COMB.pickle 
                          --fasta sample_10M_genes.fastq.gz 
                          --type_factorization ICFL_COMB -n 4

          - `return`    : creates the file test_rf_fingerprint_result.txt containing 
                          a row for each read in the FASTA file. 
        ```
        
    `Rule-based read classifier:`
    
        ```
        - `method`    : testing_reads_majority_step in SCRIPT testing.py

        - `cmd line`  : python testing.py --step test_majority --path testing/ 
                        --fasta sample_10M_genes.fastq.gz --best_model RF_ICFL_COMB_K5.pickle 
                        --criterion majority --type_factorization ICFL_COMB --k_value 5 -n 4

        - `return`    : creates a file test_majority_result.txt containing a row 
                        for each read in the FASTA file. 
        ```

- `Compute metrics:`

       ``` 
       - `cmd line`   : python metrics.py --path fingerprint/test/ 
                        --file test_majority_result_no_thresholds_list.txt 
                        --problem classification
       ```
