# GAME-consistency-evaluators

To evaluate whether models were internally consistent in their predictions, we created Consistency Evaluators. We randomly sampled 150 ATAC-seq peaks from an iPSC cell line (ENCFF121CAA) and 100kbp sequences at the center of each ATAC-seq peak were extracted. Each of the 150 sequences were mono- and dinucleotide- shuffled using Biasaway, and each of their corresponding reverse complements were included in the Evaluator data for a total of 900 sequences. The final sequences used in these Evaluators are in `evaluator_data/all_consistency_data.csv`

Using these same set of sequences we designed three Evaluators, two for K562 cell type in Homo sapiens (point and track request) and the third for "Astro" (Astrocytes) cell type in Mus musculus (point request). Users can re-build similar Evaluators for any cell type of interes using the code in the corresponding folders. For point requests, the Pearson correlation (r) was calculated between the single predicted values for forward and reverse complement sequence predictions. For track requests, the Evaluator calculates the Pearson correlation between the forward predictions and the reversed reverse complement predictions for each sequence. The final metric is reported as the mean of the correlations across pairs. 

## Usage: 

Prebuilt containers can be downloaded from Hugging Face:

- Consistency Point Request (K562, homo_sapiens): https://huggingface.co/datasets/deBoerLab/Consistency_Point_K562_GAME
- Consistency Track Request (K562, homo_sapiens): https://huggingface.co/datasets/deBoerLab/Consistency_Track__GAME
- Consistency Point Request (Astro, mus_musculus): https://huggingface.co/datasets/deBoerLab/Consistency_Point_Astro_GAME

```bash
apptainer run --containall \
    -B /path/to/evaluator_data:/evaluator_data \
    -B /path/to/predictions:/predictions \
    consistency_evaluator_{}.sif PREDICTOR_HOST PREDICTOR_PORT /predictions
```

## Development
These codebases serve as a templates for creating custom evaluators that follow a similar structure. The modular design makes it easy to adapt for different evaluation tasks.